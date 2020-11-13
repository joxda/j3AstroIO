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

#include "fitsio.h"
#include "libraw/libraw.h"
#include <exiv2/exiv2.hpp>

#include <magic.h>
#include <iostream>
#include <stdio.h>

#include "opencv2/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include "j3AstroIO.hpp"

typedef Exiv2::ExifData::const_iterator (*EasyAccessFct)(
    const Exiv2::ExifData &ed);

struct EasyAccess
{
    const char* label_;
    EasyAccessFct findFct_;
};

void printerror(int status)
{
    /* Print out cfitsio error messages and exit program */
    if (status)
    {
        fits_report_error(stderr, status); /* print error report */

        exit(status); /* terminate the program, returning error status */
    }
    return;
}

int writeTif(const char* ofile, cv::InputArray output, float factor)
{
    // TBD 8bit 16bit 32bit? jpg png?? TIF COMPRESSION? EXIF DATA?
    cv::Mat out;
    if (output.channels() == 3)
    {
        output.getMat().convertTo(out, CV_16UC3, factor);
    }
    if (output.channels() == 1)
    {
        output.getMat().convertTo(out, CV_16UC1, factor);
    }

    cv::imwrite(ofile, out);
    return 0;
}

int writeFits(const char* ofile, cv::InputArray outA)
{
    cv::Mat out = outA.getMat();
    //cv::Mat outp(cv::Size(out.size().width, out.size().height), CV_32FC1);
    cv::Mat output; //(cv::Size(out.size().height, out.size().width), out.type());
    //out.convertTo(outp, CV_32FC1);
    flip(out, output, 0); // TBD doe away with flip by indexing properly

    fitsfile* fptr;
    int status, ii, jj, bitpix;
    long fpixel, nelements;

    int naxis = 2;
    long naxes[2] = {output.size().width,
                     output.size().height
                    };

    float* farray[output.size().height];
    double* darray[output.size().height];
    unsigned short* usarray[output.size().height];
    short* sarray[output.size().height];
    int* iarray[output.size().height];
    unsigned char* ucarray[output.size().height];
    char* carray[output.size().height];

    std::cout << out.channels() << std::endl;

    remove(ofile); /* Delete old file if it already exists */

    status = 0; /* initialize status before calling fitsio routines */

    if (fits_create_file(&fptr, ofile, &status)) /* create new FITS file */
        printerror(status); /* call printerror if error occurs */

    /* write the required keywords for the primary array image.     */
    /* Since bitpix = USHORT_IMG, this will cause cfitsio to create */
    /* a FITS image with BITPIX = 16 (signed short integers) with   */
    /* BSCALE = 1.0 and BZERO = 32768.  This is the convention that */
    /* FITS uses to store unsigned integers.  Note that the BSCALE  */
    /* and BZERO keywords will be automatically written by cfitsio  */
    /* in this case.                                                */

    //int na3 = 3;
    std::vector<cv::Mat> planes;
    if(out.channels() > 1)
    {
        if (fits_create_img(fptr, 8, 0, 0, &status))
            printerror(status);
        cv::split(out, planes);
    }
    else
    {
        planes.push_back(out);
    }
    for(int c = planes.size() - 1; c >= 0; c--)
    {

        switch (out.type() & CV_MAT_DEPTH_MASK)
        {
            case CV_64F:
                std::cout << "CV_64F" << std::endl;
                if (fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status))
                    printerror(status);
                darray[0] = (double*)malloc(naxes[0] * naxes[1] * sizeof(double));
                /* initialize pointers to the start of each row of the image */
                for (ii = 1; ii < naxes[1]; ii++)
                    darray[ii] = darray[ii - 1] + naxes[0];
                for (ii = 0; ii < naxes[0]; ii++)
                {
                    for (jj = 0; jj < naxes[1]; jj++)
                    {
                        darray[jj][ii] = *planes[c].ptr<double>(jj, ii);
                    }
                }
                fpixel = 1;                      /* first pixel to write      */
                nelements = naxes[0] * naxes[1]; /* number of pixels to write */

                /* write the array of unsigned integers to the FITS file */
                if (fits_write_img(fptr, TDOUBLE, fpixel, nelements, darray[0], &status))
                    printerror(status);
                //free(darray[0]); /* free previously allocated memory */
                break;
            case CV_32F:
                std::cout << "CV_32F" << std::endl;
                if (fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status))
                    printerror(status);
                farray[0] = (float*)malloc(naxes[0] * naxes[1] * sizeof(float));
                /* initialize pointers to the start of each row of the image */
                for (ii = 1; ii < naxes[1]; ii++)
                    farray[ii] = farray[ii - 1] + naxes[0];
                for (ii = 0; ii < naxes[0]; ii++)
                {
                    for (jj = 0; jj < naxes[1]; jj++)
                    {
                        farray[jj][ii] = *planes[c].ptr<float>(jj, ii);
                    }
                }
                fpixel = 1;                      /* first pixel to write      */
                nelements = naxes[0] * naxes[1]; /* number of pixels to write */

                /* write the array of unsigned integers to the FITS file */
                if (fits_write_img(fptr, TFLOAT, fpixel, nelements, farray[0], &status))
                    printerror(status);
                //free(farray[0]); /* free previously allocated memory */
                break;
            case CV_32S:
                std::cout << "CV_32S" << std::endl;
                if (fits_create_img(fptr, LONG_IMG, naxis, naxes, &status))
                    printerror(status);
                iarray[0] = (int*)malloc(naxes[0] * naxes[1] * sizeof(int));
                /* initialize pointers to the start of each row of the image */
                for (ii = 1; ii < naxes[1]; ii++)
                    iarray[ii] = iarray[ii - 1] + naxes[0];

                for (ii = 0; ii < naxes[0]; ii++)
                {
                    for (jj = 0; jj < naxes[1]; jj++)
                    {
                        iarray[jj][ii] = *planes[c].ptr<int>(jj, ii);
                    }
                }
                fpixel = 1;                      /* first pixel to write      */
                nelements = naxes[0] * naxes[1]; /* number of pixels to write */

                /* write the array of unsigned integers to the FITS file */
                if (fits_write_img(fptr, TLONG, fpixel, nelements, iarray[0], &status))
                    printerror(status);
                //free(iarray[0]); /* free previously allocated memory */
                break;
            case CV_16U:
                std::cout << "CV_16U" << std::endl;
                if (fits_create_img(fptr, USHORT_IMG, naxis, naxes, &status))
                    printerror(status);
                usarray[0] = (unsigned short*)malloc(naxes[0] * naxes[1] * sizeof(unsigned short));
                /* initialize pointers to the start of each row of the image */
                for (ii = 1; ii < naxes[1]; ii++)
                    usarray[ii] = usarray[ii - 1] + naxes[0];
                for (ii = 0; ii < naxes[0]; ii++)
                {
                    for (jj = 0; jj < naxes[1]; jj++)
                    {
                        usarray[jj][ii] = *planes[c].ptr<unsigned short>(jj, ii);
                    }
                }
                fpixel = 1;                      /* first pixel to write      */
                nelements = naxes[0] * naxes[1]; /* number of pixels to write */

                /* write the array of unsigned integers to the FITS file */

                std::cout << "T" << std::endl;
                if (fits_write_img(fptr, TUSHORT, fpixel, nelements, usarray[0], &status))
                    printerror(status);
                std::cout << "T" << std::endl;
                break;
            case CV_16S:
                std::cout << "CV_16S" << std::endl;
                if (fits_create_img(fptr, SHORT_IMG, naxis, naxes, &status))
                    printerror(status);
                sarray[0] = (short*)malloc(naxes[0] * naxes[1] * sizeof(short));
                /* initialize pointers to the start of each row of the image */
                for (ii = 1; ii < naxes[1]; ii++)
                    sarray[ii] = sarray[ii - 1] + naxes[0];
                for (ii = 0; ii < naxes[0]; ii++)
                {
                    for (jj = 0; jj < naxes[1]; jj++)
                    {
                        sarray[jj][ii] = *planes[c].ptr<short>(jj, ii);
                    }
                }
                fpixel = 1;                      /* first pixel to write      */
                nelements = naxes[0] * naxes[1]; /* number of pixels to write */

                /* write the array of unsigned integers to the FITS file */
                if (fits_write_img(fptr, TSHORT, fpixel, nelements, sarray[0], &status))
                    printerror(status);
                break;
            case CV_8S:
                std::cout << "CV_8S" << std::endl;
                if (fits_create_img(fptr, SBYTE_IMG, naxis, naxes, &status))
                    printerror(status);
                carray[0] = (char*)malloc(naxes[0] * naxes[1] * sizeof(char));
                /* initialize pointers to the start of each row of the image */
                for (ii = 1; ii < naxes[1]; ii++)
                    carray[ii] = carray[ii - 1] + naxes[0];

                for (ii = 0; ii < naxes[0]; ii++)
                {
                    for (jj = 0; jj < naxes[1]; jj++)
                    {
                        carray[jj][ii] = *planes[c].ptr<char>(jj, ii);
                    }
                }
                fpixel = 1;                      /* first pixel to write      */
                nelements = naxes[0] * naxes[1]; /* number of pixels to write */

                /* write the array of unsigned integers to the FITS file */
                if (fits_write_img(fptr, TSBYTE, fpixel, nelements, carray[0], &status))
                    printerror(status);
                break;
            case CV_8U:

                std::cout << "CV_8U" << std::endl;
                if (fits_create_img(fptr, BYTE_IMG, naxis, naxes, &status))
                    printerror(status);
                ucarray[0] = (unsigned char*)malloc(naxes[0] * naxes[1] * sizeof(unsigned char));
                /* initialize pointers to the start of each row of the image */
                for (ii = 1; ii < naxes[1]; ii++)
                    ucarray[ii] = ucarray[ii - 1] + naxes[0];
                for (ii = 0; ii < naxes[0]; ii++)
                {
                    for (jj = 0; jj < naxes[1]; jj++)
                    {
                        ucarray[jj][ii] = *planes[c].ptr<unsigned char>(jj, ii);
                    }
                }
                fpixel = 1;                      /* first pixel to write      */
                nelements = naxes[0] * naxes[1]; /* number of pixels to write */

                /* write the array of unsigned integers to the FITS file */
                if (fits_write_img(fptr, TBYTE, fpixel, nelements, ucarray[0], &status))
                    printerror(status);
                break;
            default:
                return 1;
        }
    }

    /* write another optional keyword to the header */
    /* Note that the ADDRESS of the value is passed in the routine */

    /*fits_open_file(&fptr, "WFPC2ASSNu5780205bx.fits", READONLY, &status);
    long naxes[2];
    fits_get_img_size(fptr, 3, naxes, &status);

    fitsfile *ofptr;
    fits_create_file(&ofptr, "o_nasa.fits", &status);
    fits_copy_header(fptr, ofptr, &status);
    */

    //  exposure = 1500;

    //            fptr, TLONG, "EXPOSURE", &exposure, "Total Exposure Time", &status))
    printerror(status);

    if (fits_close_file(fptr, &status)) /* close the file */
        printerror(status);

    return 0;
}


int writeFile(char* ofile, cv::InputArray output)
{
    char* ext;
    ext = std::strrchr(ofile, '.');
    if (std::strcmp(ext, ".fits") == 0 || std::strcmp(ext, ".fit") == 0)
    {
        writeFits(ofile, output);
    }
    else
        //if (strcmp(ext, ".tiff") == 0 || strcmp(ext, ".tif") == 0)
    {
        writeTif(ofile, output);
    }
    return 0;
}



int open(const char* file, cv::OutputArray image)
{

    std::string str = mime(file);

    int success = -1;
    // TBD: CHECK WHETHER ACTUALLY SUCCESSFUL
    // TBD: ERROR MESSAGE AND GRACEFUL EXIT IF NOT

    if (str.find("image/fits") == 0)
    {
        open_fits(file, image);
        success = 0;
    }
    else if (str.find("image/x") == 0)
    {
        std::cout << "raw" << std::endl;
        open_raw(file, image);
        success = 0;
    }
    else if (cv::haveImageReader(file)) // TBD!! not in v3
    {
        success = open_opencv(file, image);
    }
    return success;
}

LensPars getPars(const char* file)
{
    LensPars par;
    try
    {
        Exiv2::XmpParser::initialize();
        ::atexit(Exiv2::XmpParser::terminate);

        // Exiv2::Image::AutoPtr EXimage = Exiv2::ImageFactory::open(file);
        Exiv2::Image::UniquePtr EXimage = Exiv2::ImageFactory::open(file);
        assert(EXimage.get() != 0);
        EXimage->readMetadata();
        Exiv2::ExifData &ed = EXimage->exifData();

        // static const EasyAccess easyAccess = {"Lens name",
        // Exiv2::lensName     }; Exiv2::ExifData::const_iterator pos =
        // easyAccess.findFct_(ed);
        Exiv2::ExifData::const_iterator pos = lensName(ed);
        Exiv2::ExifData::const_iterator posFN = fNumber(ed);
        Exiv2::ExifData::const_iterator posFL = focalLength(ed);
        // lensName exposureTime focalLength
        // std::cout << std::setw(20) << std::left << easyAccess.label_;
        /*if (pos != ed.end()) {
            std::cout << " (" << std::setw(35) << pos->key() << ") : "
                          << pos->print(&ed) << "\n";
        } else {
                std::cout << " (" << std::setw(35) << " " << ") : \n";
            }
          */

        Exiv2::ExifKey key("Exif.Photo.FocalPlaneXResolution");
        Exiv2::ExifData::iterator CRpos = ed.findKey(key);
        if (CRpos == ed.end())
            std::cout << "THIS IS THE PROBLEM" << std::endl;
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

        par.name = std::string(pos->print(&ed));
        par.apertureN = std::stof((std::string(posFN->print(&ed))).substr(1, 10));
        par.focalLength = std::stof(posFL->print(&ed));
    }
    catch (Exiv2::AnyError &e)
    {
        std::cout << "Caught Exiv2 exception '" << e << "'\n";
    }
    return par;
}


int open_raw(const char* file, cv::OutputArray image)
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
        cv::Mat(cv::Size(imag->width, imag->height), CV_16UC3, imag->data, cv::Mat::AUTO_STEP)
        .clone();

    LibRaw::dcraw_clear_mem(imag);
    cvtColor(im, image, cv::COLOR_RGB2BGR);

    iProcessor.recycle();
    return ret;
}

int open_fits(const char* file, cv::OutputArray image)
{
    fitsfile* fptr;

    int status = 0; /* CFITSIO status value MUST be initialized to zero! */
    int hdutype, naxis, datatype = 0, anynul;
    long naxes[2], totpix;
    float *array, bscale = 1.0, bzero = 0.0, nulval = 0.;

    int cvtype;
    if (!fits_open_image(&fptr, file, READONLY, &status))
    {
        if (fits_get_hdu_type(fptr, &hdutype, &status) || hdutype != IMAGE_HDU)
        {
            printf("Error: this program only works on images, not tables\n");
            return 1;
        }

        fits_get_img_dim(fptr, &naxis, &status);
        fits_get_img_size(fptr, 2, naxes, &status);

        if (status || naxis != 2)
        {
            printf(
                "Error: NAXIS = %d.  Only 2-D images are supported.\n", naxis);
            return 1;
        }

        int bitpix, bytepix;

        fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status);
        switch (bitpix) // multi-extension fits files for color?
        {
            case BYTE_IMG:
                datatype = TBYTE;
                cvtype = CV_8UC1;
                // sep_type = SEP_TBYTE;
                break;
            case SBYTE_IMG:
                datatype = TBYTE;
                cvtype = CV_8SC1;
                // sep_type = SEP_TBYTE;
                break;
            case SHORT_IMG:
                datatype = TSHORT;
                cvtype = CV_16SC1;
                break;
            case USHORT_IMG:
                datatype = TSHORT;
                cvtype = CV_16UC1;
                break;
            case LONG_IMG:
                datatype = TINT;
                cvtype = CV_32SC1;
                //       sep_type = SEP_TINT ;
                break;
            case FLOAT_IMG:
                datatype = TFLOAT;
                cvtype = CV_32FC1;
                //           sep_type = SEP_TFLOAT ;
                break;
            case DOUBLE_IMG:
                datatype = TDOUBLE;
                cvtype = CV_64FC1;
                //           sep_type = SEP_TDOUBLE;
                break;
            default:
                return 1;
        }

        bytepix = abs(bitpix) / 8;
        totpix = naxes[0] * naxes[1];
        int width = (int)naxes[0];
        int height = (int)naxes[1];

        array = (float*)calloc(totpix, bytepix); // TBD will this work also for eg 8bit data? double?
        if (!array)
        {
            printf("Memory allocation error\n");
            return (0);
        }

        fits_set_bscale(fptr, bscale, bzero, &status);

        long first = 1;
        //        while (totpix > 0 && !status)
        //       {
        /* read all or part of image then write it back to the output file */
        /// SEE ORIGINAL SAMPLE CODE FOR READING IN CHUNKS?
        fits_read_img(
            fptr, datatype, first, totpix, &nulval, array, &anynul, &status);

        fits_close_file(fptr, &status);

        cv::Mat fitin =
            cv::Mat(cv::Size(width, height), cvtype, array, cv::Mat::AUTO_STEP);

        flip(fitin, image, 0);
    }
    printf("DONE FITS\n");
    return status;
}


int open_opencv(const char* file, cv::OutputArray image)
{
    cv::Mat im = image.getMat();
    im = cv::imread(file, cv::IMREAD_COLOR | cv::IMREAD_ANYDEPTH);
    if (image.empty())
        return -1;
    // Read the file
    // open some other formates with opencv
    return 0;
}

std::string mime(const char* file)
{
    struct magic_set* myt = magic_open(MAGIC_CONTINUE |
                                       MAGIC_ERROR /*|MAGIC_DEBUG*/ | MAGIC_PRESERVE_ATIME | MAGIC_MIME);
    // magic_t myt =
    // magic_open(MAGIC_CONTINUE|MAGIC_ERROR/*|MAGIC_DEBUG*/|MAGIC_MIME_TYPE);
    if (myt == NULL)
    {
        printf("ERROR\n");
    }
    // printf("%s\n",magic_version());

    int status = magic_load(myt,
                            NULL /*"./magic.mgc"*/); // TBD do this copy thing and get the path
    // relative to project...
    // TBD if not == 0 -> error with magic.mgc...
    if (status != 0)
    {
        printf("ERROR\n");
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




