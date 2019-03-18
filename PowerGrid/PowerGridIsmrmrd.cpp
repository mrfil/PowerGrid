/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [PowerGridIsmrmrd.cpp]

    Synopsis    [PowerGrid reconstruction executable supporting ISMRMRD format
                                as input.]

    Description [This reconstruction supports 2D and 3D reconstructions with
                                time segmentation for field correction. This
 support is experimental.]

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

// //Project headers.
#include "PowerGrid.h"
#include "processIsmrmrd.hpp"
#include "processNIFTI.hpp"
#include <boost/program_options.hpp>
#include <chrono>
namespace po = boost::program_options;
//using namespace PowerGrid;

int main(int argc, char **argv) {
  std::string rawDataFilePath, outputImageFilePath, senseMapFilePath,
      fieldMapFilePath, precisionString, TimeSegmentationInterp, FourierTrans,
      rawDataNavFilePath;
  uword Nx, Ny, Nz, NShots = 1, type = 1, L = 0, NIter = 10, FtType = 0;
  //uword ;
  double beta = 0.0;
  uword dims2penalize = 3;
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "produce help message")(
      "inputData,i", po::value<std::string>(&rawDataFilePath)->required(),
      "input ISMRMRD Raw Data file")
 			("outputImage,o", po::value<std::string>(&outputImageFilePath)->required(), "output file path for NIFTIimages")

      /*
      ("inputDataNav,-N", po::value<std::string>(&rawDataNavFilePath), "input
      ISMRMRD Navigator Raw Data")
      ("outputImage,o",
      po::value<std::string>(&outputImageFilePath)->required(), "output ISMRMRD
      Image file")
      ("SENSEMap,S", po::value<std::string>(&senseMapFilePath),
       "Enable SENSE recon with the specified SENSE map in ISMRMRD image
      format")
      ("FieldMap,F", po::value<std::string>(&fieldMapFilePath),
      "Enable field corrected reconstruction with the specified field map in
      ISMRMRD format")
      ("Precision,P", po::value<std::string>(&precisionString),
       "Numerical precision to use, float or double currently supported")
       */
      ("Nx,x", po::value<uword>(&Nx)->required(), "Image size in X (Required)")
          ("Ny,y", po::value<uword>(&Ny)->required(), "Image size in Y (Required)")
          ("Nz,z", po::value<uword>(&Nz)->required(), "Image size in Z (Required)")
          ("NShots,s", po::value<uword>(&NShots), "Number of shots per image")
          ("TimeSegmentationInterp,I", po::value<std::string>(&TimeSegmentationInterp)->required(), "Field Correction Interpolator (Required)")
          ("FourierTransform,F", po::value<std::string>(&FourierTrans)->required(), "Implementation of Fourier Transform")
          ("TimeSegments,t", po::value<uword>(&L)->required(), "Number of time segments (Required)")
          ("Beta,B", po::value<double>(&beta), "Spatial regularization penalty weight")
          ("CGIterations,n", po::value<uword>(&NIter), "Number of preconditioned conjugate gradient interations for main solver")
          ("Dims2Penalize,D", po::value<uword>(&dims2penalize), "Dimensions to apply regularization to (2 or 3).");


  po::variables_map vm;

  try {

    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cout << desc << std::endl;
      return 1;
    }
    /*
    if(precisionString.compare("double") ==0) {
            typedef double PGPrecision;
    } else if(precisionString.compare("float") == 0) {
            typedef float PGPrecision;
    } else {
            typedef double PGPrecision;
            std::cout << "Did not recognize precision option. Defaulting to
    double precision." << std::endl;
    }
    */

    if (FourierTrans.compare("DFT") == 0) {
        FtType = 2;
    
    } else if (FourierTrans.compare("NUFFT") == 0) {
      
      FtType = 1;
      if (TimeSegmentationInterp.compare("hanning") == 0) {
        type = 1;
      } else if (TimeSegmentationInterp.compare("minmax") == 0) {
        type = 2;
      } else if (TimeSegmentationInterp.compare("histo") == 0) {
        type = 3;
      } else {
        std::cout << "Did not recognize temporal interpolator selection. " << std::endl
                  << "Acceptable values are hanning or minmax."            << std::endl;
        return 1;
      }
    } else if (FourierTrans.compare("DFTGrads") == 0) {
      FtType = 3;
    } else {
      std::cout << "Did not recognize Fourier transform selection. " << std::endl
                << "Acceptable values are DFT or NUFFT."             << std::endl;
      return 1;
    }


  } catch (boost::program_options::error &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }

  Col<float> FM;
  Col<float> fmSlice;
  Col<std::complex<float>> sen;
 	Col<std::complex<float>> senSlice;

  ISMRMRD::Dataset *d;
  ISMRMRD::IsmrmrdHeader hdr;
	acqTracking *acqTrack;
  processISMRMRDInput<float>(rawDataFilePath, d, hdr, FM, sen, acqTrack);

  std::cout << "Number of elements in SENSE Map = " << sen.n_rows << std::endl;
  std::cout << "Number of elements in Field Map = " << FM.n_rows << std::endl;

  uword numAcq = d->getNumberOfAcquisitions();

  // Grab first acquisition to get parameters (We assume all subsequent
  // acquisitions will be similar).
  ISMRMRD::Acquisition acq;
  d->readAcquisition(0, acq);
  uword nro = acq.number_of_samples();
  uword nc = acq.active_channels();

  Col<float> ix, iy, iz;
  initImageSpaceCoords(ix, iy, iz, Nx, Ny, Nz);
  Col<float> kx(nro), ky(nro), kz(nro), tvec(nro);
  Col<std::complex<float>> data(nro * nc);
  Col<std::complex<float>> ImageTemp(Nx * Ny * Nz);

  // Check and abort if we have more than one encoding space (No Navigators for
  // now).
  std::cout << "hdr.encoding.size() = " << hdr.encoding.size() << std::endl;
  if (hdr.encoding.size() != 1) {
    std::cout << "There are " << hdr.encoding.size()
              << " encoding spaces in this file" << std::endl;
    std::cout
        << "This recon does not handle more than one encoding space. Aborting."
        << std::endl;
    return -1;
  }

  int NShotMax  = hdr.encoding[0].encodingLimits.kspace_encoding_step_1->maximum;
  int NParMax   = hdr.encoding[0].encodingLimits.kspace_encoding_step_2->maximum;
  int NSliceMax = hdr.encoding[0].encodingLimits.slice->maximum;
  int NSetMax   = hdr.encoding[0].encodingLimits.set->maximum;
  int NRepMax   = hdr.encoding[0].encodingLimits.repetition->maximum;
  int NAvgMax   = hdr.encoding[0].encodingLimits.average->maximum;
  int NSegMax   = hdr.encoding[0].encodingLimits.segment->maximum;
  int NEchoMax  = hdr.encoding[0].encodingLimits.contrast->maximum;
  int NPhaseMax = hdr.encoding[0].encodingLimits.phase->maximum;

  std::cout << "NParMax = "   << NParMax << std::endl;
  std::cout << "NShotMax = "  << NShotMax << std::endl;
  std::cout << "NSliceMax = " << NSliceMax << std::endl;
  std::cout << "NSetMax = "   << NSetMax << std::endl;
  std::cout << "NRepMax = "   << NRepMax << std::endl;
  std::cout << "NAvgMax = "   << NAvgMax << std::endl;
  std::cout << "NEchoMax = "  << NEchoMax << std::endl;
  std::cout << "NPhaseMax = " << NPhaseMax << std::endl;
  std::cout << "NSegMax = "   << NSegMax << std::endl;

   std::cout << "About to loop through the counters and scan the file"
            << std::endl;

  std::string baseFilename = "img";
  std::string filename;
  if (!outputImageFilePath.empty() && *outputImageFilePath.rbegin() != '/') {
    	outputImageFilePath += '/';
	}
  //uword NSet = 0; //Set is only used for arrayed ADCs
  //uword NSeg = 0;
	for (uword NPhase = 0; NPhase <= NPhaseMax; NPhase++) {
		for (uword NEcho = 0; NEcho <= NEchoMax; NEcho++) {
            for (uword NAvg = 0; NAvg <= NAvgMax; NAvg++) {
                for (uword NRep = 0; NRep < NRepMax +1; NRep++) {
                    for (uword NSlice = 0; NSlice<=NSliceMax; NSlice++) {
                      filename = outputImageFilePath + baseFilename + "_" + "Slice" + std::to_string(NSlice) +
          								"_" + "Rep" + std::to_string(NRep) + "_" + "Avg" + std::to_string(NAvg) +
          								"_" + "Echo" + std::to_string(NEcho) + "_" + "Phase" + std::to_string(NPhase);
                          
	                    senSlice = getISMRMRDCompleteSENSEMap<std::complex<float>>(d, sen, NSlice, Nx*Ny*Nz);
						          fmSlice = getISMRMRDCompleteFieldMap<float>(d, FM, NSlice, (uword) (Nx*Ny*Nz));

	                    getCompleteISMRMRDAcqData<float>(d, acqTrack, NSlice, NRep, NAvg, NEcho, NPhase, data, kx, ky,
			                    kz, tvec);
	                    std::cout << "Number of elements in kx = " << kx.n_rows << std::endl;
	                    std::cout << "Number of elements in ky = " << ky.n_rows << std::endl;
	                    std::cout << "Number of elements in kz = " << kz.n_rows << std::endl;
	                    std::cout << "Number of rows in data = " << data.n_rows << std::endl;
	                    std::cout << "Number of columns in data = " << data.n_cols << std::endl;

	                    QuadPenalty<float> R(Nx, Ny, Nz, beta, dims2penalize);

                      if (FtType == 1) {
	                      Gnufft<float> G(kx.n_rows, (float) 2.0, Nx, Ny, Nz, kx, ky, kz, ix,
			                    iy, iz);
	                      TimeSegmentation<float, Gnufft<float>> A(G, fmSlice, tvec, kx.n_rows, Nx*Ny*Nz, L, type, NShots);
                        SENSE<float, TimeSegmentation<float, Gnufft<float>>> Sg(A, senSlice, kx.n_rows, Nx*Ny*Nz, nc);
	                      ImageTemp = reconSolve<float, SENSE<float, TimeSegmentation<float, Gnufft<float>>>,
			                    QuadPenalty<float>>(data, Sg, R, kx, ky, kz, Nx,
			                    Ny, Nz, tvec, NIter);
                      } else if (FtType == 2) {
                        Gdft<float> A(kx.n_rows, Nx*Ny*Nz,kx,ky,kz,ix,iy,iz,fmSlice,tvec);
	                      SENSE<float, Gdft<float>> Sg(A, senSlice, kx.n_rows, Nx*Ny*Nz, nc);
	                      ImageTemp = reconSolve<float, SENSE<float, Gdft<float>>,
	                            QuadPenalty<float>>(data, Sg, R, kx, ky, kz, Nx,
                              Ny, Nz, tvec, NIter);
                      } else if (FtType == 3) {
                        GdftR2<float> A(kx.n_rows, Nx*Ny*Nz,kx,ky,kz,ix,iy,iz,fmSlice,tvec,Nx,Ny,Nz);
	                      SENSE<float, GdftR2<float>> Sg(A, senSlice, kx.n_rows, Nx*Ny*Nz, nc);
	                      ImageTemp = reconSolve<float, SENSE<float, GdftR2<float>>,
	                            QuadPenalty<float>>(data, Sg, R, kx, ky, kz, Nx,
                              Ny, Nz, tvec, NIter);
                      }


	                  //writeISMRMRDImageData<float>(d, ImageTemp, Nx, Ny, Nz);
                    writeNiftiMagPhsImage<float>(filename,ImageTemp,Nx,Ny,Nz);
                    }

                }
            }
        }
    }


  // Close ISMRMRD::Dataset, hdr, and acqTrack
	closeISMRMRDData(d,hdr,acqTrack);

  return 0;
}
