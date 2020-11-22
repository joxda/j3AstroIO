/*******************************************************************************
  Copyright(c) 2020 Joachim Janz. All rights reserved.

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the Free
  Software Foundation; either version 2 of the License, or (at your option)
  any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU Library General Public License
  along with this library; see the file COPYING.LIB.  If not, write to
  the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
  Boston, MA 02110-1301, USA.

  The full GNU General Public License is included in this distribution in the
  file called LICENSE.

*******************************************************************************/
//
//  j3AstroIO.cpp
//
//  Created by Joachim Janz on 12/11/2020.
//  Copyright © 2020 Joachim Janz. All rights reserved.
//

#include <exiv2/exiv2.hpp>
#include <exiv2/exv_conf.h>
#ifdef EXIV2_VERSION
# ifndef EXIV2_TEST_VERSION
# define EXIV2_TEST_VERSION(major,minor,patch) \
    ( EXIV2_VERSION >= EXIV2_MAKE_VERSION(major,minor,patch) )
# endif
#else
# define EXIV2_TEST_VERSION(major,minor,patch) (false)
#endif

#include <iostream>
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include "j3AstroIOexiv.hpp"

typedef Exiv2::ExifData::const_iterator (*EasyAccessFct)(
    const Exiv2::ExifData &ed);

struct EasyAccess
{
    const char* label_;
    EasyAccessFct findFct_;
};



PhotoPars getPars(const char* file)
{
    PhotoPars par;
    try
    {
        Exiv2::XmpParser::initialize();
        ::atexit(Exiv2::XmpParser::terminate);
        
        #if EXIV2_TEST_VERSION(0, 27, 1)
            Exiv2::Image::AutoPtr EXimage = Exiv2::ImageFactory::open(file);
        #else
            Exiv2::Image::UniquePtr EXimage = Exiv2::ImageFactory::open(file);
        #endif
        assert(EXimage.get() != 0);
        EXimage->readMetadata();
        Exiv2::ExifData &ed = EXimage->exifData();

        Exiv2::ExifData::const_iterator pos = lensName(ed);
        Exiv2::ExifData::const_iterator posFN = fNumber(ed);
        Exiv2::ExifData::const_iterator posFL = focalLength(ed);
        Exiv2::ExifData::const_iterator posCamMake = make(ed);
        Exiv2::ExifData::const_iterator posCamModel = model(ed);

        Exiv2::ExifKey key("Exif.Photo.FocalPlaneXResolution");
        Exiv2::ExifData::iterator CRpos = ed.findKey(key);
        if (CRpos == ed.end())
            std::cout << "Did not find key" << std::endl;
        float xres = std::stof(CRpos->print(&ed));
        key = Exiv2::ExifKey("Exif.Photo.FocalPlaneYResolution");
        CRpos = ed.findKey(key);
        float yres = std::stof(CRpos->print(&ed));
        key = Exiv2::ExifKey("Exif.Photo.FocalPlaneResolutionUnit");
        CRpos = ed.findKey(key);
        float fac = CRpos->print(&ed) == "inch" ? 25.4 : 10.;
        xres /= fac;
        yres /= fac;
        float xmm, ymm;
        key = Exiv2::ExifKey("Exif.Photo.PixelXDimension");
        CRpos = ed.findKey(key);
        xmm = std::stof(CRpos->print(&ed)) / xres;
        key = Exiv2::ExifKey("Exif.Photo.PixelYDimension");
        CRpos = ed.findKey(key);
        ymm = std::stof(CRpos->print(&ed)) / yres;
        par.cropFactor = sqrt(xmm * xmm + ymm * ymm) / 43.267;

        par.lensName = std::string(pos->print(&ed));
        par.camName = std::string(posCamModelMake->print(&ed));
        par.camMake = std::string(posCamMake->print(&ed));
        par.apertureN = std::stof((std::string(posFN->print(&ed))).substr(1, 10));
        par.focalLength = std::stof(posFL->print(&ed));
    }
    catch (Exiv2::AnyError &e)
    {
        std::cout << "Caught Exiv2 exception '" << e << "'\n";
    }
    return par;
}


