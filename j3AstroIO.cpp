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

#include <magic.h>

#ifdef CIRUN
  #include "resource.h"
#endif

#include <cstring>

#include <iostream>
#include <stdio.h>
#include <strings.h>


#include "j3AstroIO.hpp"


int open(const char* file)
{
    std::string str = mime(file);
    std::cout << "MIME: " << file << " " << str.c_str() << std::endl;
    int success = 0;

    return success;
}

std::string mime(const char* file)
{
    struct magic_set* myt = magic_open(MAGIC_CONTINUE |
                                       MAGIC_ERROR /*|MAGIC_DEBUG*/ | MAGIC_PRESERVE_ATIME | MAGIC_MIME);
    // magic_t myt =
    // magic_open(MAGIC_CONTINUE|MAGIC_ERROR/*|MAGIC_DEBUG*/|MAGIC_MIME_TYPE);
    if (myt == NULL)
    {
        printf("Magic open ERROR\n");
    }
    // printf("%s\n",magic_version());
#ifdef CIRUN
    #pragma message("USE EMBEDDED MAGIC")
    Resource text = LOAD_RESOURCE(magic_mgc);
    std::cout << "LOADING BINARY MAGIC" << std::endl;
    int status = magic_load_buffers(myt, (void**)&text.data(), (size_t*)&text.size(), 1);
    std::cout << status << std::endl;
#else
    int status = magic_load(myt,
                            NULL /*"./magic.mgc"*/); // TBD do this copy thing and get the path
#endif
    // relative to project...
    // TBD if not == 0 -> error with magic.mgc...
    if (status != 0)
    {
        int status = magic_load(myt,
                            "/usr/share/file/magic.mgc"/*"./magic.mgc"*/);
        while (status != 0) {
            std::string mgcfile;
            std::cout << "Magic load ERROR -please give the path to magic.mgc file: " << std::flush;
            std::cin >> mgcfile;
            status = magic_load(myt,
                            mgcfile.c_str());
            std::cout << status << std::endl;
        }
    }

    const char* mm = magic_file(myt, file);
    if (mm != NULL)
    {
        std::cout << mm << std::endl;
    }

    const unsigned long len = std::strlen(mm);

    if (mm == NULL)
    {
        (void)printf("ERROR: %s\n", magic_error(myt));
    }
    char* mimeret = new char[len + 1];

    std::strncpy(mimeret, mm, len);

    magic_close(myt);
    std::string mimereturn(mimeret);

    delete[] mimeret;
    return mimereturn;
}








