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


#include "j3AstroIOexiv.hpp"
#include <iostream>

int main(int argc, char** argv)
{
    if(argc != 2)
    {
        std::cout << "  Usage: j3AstroIOconvert inputFile" << std::endl;
        return 1;
    }
    LensPars p = getPars(argv[1]);
    std::cout << p.name << "  focal length: " << p.focalLength << " f: " << p.apertureN << " crop: " << p.cropFactor << std::endl;

    return 0;
}