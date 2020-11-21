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
//  j3AstroIOraw.cpp
//
//  Created by Joachim Janz on 12/11/2020.
//  Copyright © 2020 Joachim Janz. All rights reserved.
//

#include "libraw/libraw.h"

#include <magic.h>
#include <iostream>
#include <stdio.h>
#include <strings.h>

#include "opencv2/core/version.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include "j3AstroIOraw.hpp"


int write_opencv(const char* ofile, cv::InputArray output, float factor, int depth)
{
    int success = 0;
    if(CV_MAJOR_VERSION > 3)
    {
        if(cv::haveImageWriter(ofile))      // TBD opencv_v3 compatibility?
        {
            cv::Mat out;
            output.getMat().convertTo(out, depth, factor);
            success = cv::imwrite(ofile, out);
        }
        else
        {
            success = 1;
            const char* ext = std::strrchr(ofile, '.');;
            std::cout << "OpenCV has no writer for " << ext << " files." << std::endl;
        }
    }
    else
    {
        cv::Mat out;
        output.getMat().convertTo(out, depth, factor);
        success = cv::imwrite(ofile, out);
    }
    return success;
}


int writeFile(char* ofile, cv::InputArray output)
{
    char* ext;
    int success;
    ext = std::strrchr(ofile, '.');
    std::cout << "EXT: " << ext << std::endl;
    success = write_opencv(ofile, output);
}

int open(const char* file, cv::Mat &image)
{
    std::string str = mime(file);

    int success = -1;

    if (str.find("image/x") == 0)
    {
        std::cout << "raw" << std::endl;
        success = open_raw(file, image);
    }
    return success;
}


int open_raw(const char* file, cv::Mat &image)
{
    LibRaw iProcessor;
    int ret;
    iProcessor.imgdata.params.use_camera_wb = 1;
    iProcessor.imgdata.params.output_tiff = 1;
    iProcessor.imgdata.params.output_bps = 16; // WHY NOT -4?
    iProcessor.imgdata.params.no_auto_bright = 1;
    iProcessor.imgdata.params.output_color = 3; // 0 = raw
    iProcessor.imgdata.params.user_qual = 3;    // 3 = AHD
    iProcessor.imgdata.params.gamm[0] = 1.0;
    iProcessor.imgdata.params.gamm[1] = 1.0;

    if ((ret = iProcessor.open_file(file)) != LIBRAW_SUCCESS)
    {
        fprintf(
            stderr, "Cannot open_file %s: %s\n", file, libraw_strerror(ret));
    };

    if ((ret = iProcessor.unpack()) != LIBRAW_SUCCESS)
    {
        fprintf(stderr, "Cannot unpack %s: %s\n", file, libraw_strerror(ret));
    }

    ret = iProcessor.dcraw_process();
    if (LIBRAW_SUCCESS != ret)
    {
        fprintf(stderr, "Cannot do postpocessing on %s: %s\n", file,
                libraw_strerror(ret));
    }

    libraw_processed_image_t* imag = iProcessor.dcraw_make_mem_image(&ret);

    cv::Mat im =
        cv::Mat(cv::Size(imag->width, imag->height), CV_16UC3, imag->data, cv::Mat::AUTO_STEP);
//        .clone();

    LibRaw::dcraw_clear_mem(imag);
    cvtColor(im, image, cv::COLOR_RGB2BGR);

    iProcessor.recycle();
    return ret;
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

     int status = magic_load(myt,
                            NULL /*"./magic.mgc"*/); // TBD do this copy thing and get the path
    // relative to project...
    // TBD if not == 0 -> error with magic.mgc...
    if (status != 0)
    {
        printf("Magic load ERROR\n");
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




