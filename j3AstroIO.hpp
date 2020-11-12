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
//  j3AstroIO.hpp
//
//  Created by Joachim Janz on 12/11/2020.
//  Copyright © 2020 Joachim Janz. All rights reserved.
//

#ifndef j3AstroIO_hpp
#define j3AstroIO_hpp

#include "opencv2/core.hpp"
#include <stdio.h>

struct LensPars  
{  
    std::string name;  
    float focalLength;  
    float apertureN;  
    float cropFactor;
}; 
    void printerror(int status);

    int open(const char* file, cv::OutputArray output);
    
    int open_raw(const char* file, cv::OutputArray output);
    int open_fits(const char* file, cv::OutputArray output);
    int open_opencv(const char* file, cv::OutputArray output);


    int writeFile(char* ofile, cv::InputArray output);

    int writeTif(const char* ofile, cv::InputArray output, float factor=65355.);
    int writeFits(const char* ofile, cv::InputArray output);
    
    std::string mime(const char* file);
    LensPars getPars(const char* file);


#endif /* j3AstroIO_hpp */