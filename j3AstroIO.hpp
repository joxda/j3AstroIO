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

/** @file
 *
 * API of the j3AstroIO library
 */

#ifndef j3AstroIO_hpp
#define j3AstroIO_hpp

#include <iostream>
#include <stdio.h>
//#include "opencv2/core.hpp"
namespace cv {
class Mat;
class _InputArray;
typedef const _InputArray &InputArray;
} // namespace cv

namespace j3AstroIO {
/**
 * @brief Stucture to hold the lens and camera parameters that are read out from
 * the EXIF data
 *
 */
struct PhotoPars {
  std::string lensName;
  std::string camName;
  std::string camMake;
  float focalLength;
  float apertureN;
  float cropFactor;
};

/**
 * @brief Read the lens and camera parameters from the EXIF data
 *
 * @param[in] file File name
 * @return PhotoPars
 */
PhotoPars getPars(const char *file);

/**
 * @brief Print error code for FITS IO from cfitsio
 *
 * @param[in] status
 */
void printerror(int status);

/**
 * @brief Open image convenience function with magic to decide file format
 *
 * @param[in] file File name
 * @param[out] output cv::Mat which the image is read to
 * @return int Status (0=Success)
 */
int open(const char *file, cv::Mat &output);

/**
 * @brief Read image from FITS file
 *
 * @param[in] file File name
 * @param[out] output cv::Mat which the image is read to
 * @return int Status (0=Success)
 */
int open_fits(const char *file, cv::Mat &output);

/**
 * @brief Read image from a file with a format that can be read by OpenCV
 *
 * @param[in] file File name
 * @param[out] output cv::Mat which the image is read to
 * @return int Status (0=Success)
 */
int open_opencv(const char *file, cv::Mat &output);

/**
 * @brief Read image from file in camera raw data format
 *
 * @param[in] file File name
 * @param[out] output cv::Mat which the image is read to
 * @return int Status (0=Success)
 */
int open_raw(const char *file, cv::Mat &output);

/**
 * @brief Convenience function for writing image.
 *
 * For extensions fit/fits the cfitsio library is used to write a fits file,
 * and otherwise OpenCV handles the extension and attempts to write a file in
 * the corresponding format
 *
 * @param[in] ofile File name
 * @param[in] output cv::Mat which should be written
 * @return int Status (0=Success)
 */
int writeFile(const char *ofile, cv::InputArray output, float factor);

/**
 * @brief Write image to file using OpenCV IO
 *
 * @param[in] ofile File Name
 * @param[in] output Image to be written to the file
 * @param[in] factor The image is multiplied with this number before it is
 * written to the file
 * @param[in] depth Depth of the output image
 * @return int
 */
int write_opencv(const char *ofile, cv::InputArray output, float factor = 1.,
                 int depth = 0);

/**
 * @brief Returns the MIME type of a file (using libmagic)
 *
 * @param[in] file File Name
 * @return std::string MIME type
 */
std::string mime(const char *file);

/**
 * @brief Copy the meta data (EXIF etc) from one file to another
 *
 * @param[in] inFile File name of the input image
 * @param outFile[in]  ile name of the output image
 * @return int Status (0=Success)
 */
int copyMeta(const char *inFile, const char *outFile);

} // namespace j3AstroIO
#endif /* j3AstroIO_hpp */
