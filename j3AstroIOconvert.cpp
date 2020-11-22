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


#include "j3AstroIOraw.hpp"
#include "opencv2/core.hpp"
#include <iostream>

int main(int argc, char** argv)
{
    if(argc < 3)
    {
        std::cout << "  Usage: j3AstroIOconvert inputFile outputFile" << std::endl;
        return 1;
    }
    float factor = 1.0;
    if(argc==4) {
      factor = atof(argv[3]);
    }
    cv::Mat im;
    open(argv[1], im);
    writeFile(argv[2], im, factor);
    return 0;
}
