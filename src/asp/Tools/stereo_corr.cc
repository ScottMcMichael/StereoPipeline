// __BEGIN_LICENSE__
//  Copyright (c) 2009-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NGT platform is licensed under the Apache License, Version 2.0 (the
//  "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


/// \file stereo_corr.cc
///

#include <boost/core/null_deleter.hpp>
#include <vw/InterestPoint.h>
#include <vw/Camera/CameraTransform.h>
#include <vw/Camera/PinholeModel.h>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <vw/Stereo/CorrelationView.h>
#include <vw/Stereo/CostFunctions.h>
#include <vw/Stereo/DisparityMap.h>
#include <asp/Tools/stereo.h>
#include <asp/Core/DemDisparity.h>
#include <asp/Core/LocalHomography.h>
#include <asp/Sessions/StereoSession.h>
#include <asp/Sessions/StereoSessionPinhole.h>
#include <xercesc/util/PlatformUtils.hpp>

#include <asp/Core/InterestPointMatching.h>
#include <vw/Stereo/StereoModel.h>

#include <asp/Core/AffineEpipolar.cc>
#include <asp/Core/AffineEpipolar.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/InterestPoint/InterestData.h>

using namespace vw;
using namespace vw::stereo;
using namespace asp;
using namespace std;


//<RM>: Added DEBUG_RM flag 
#define DEBUG_RM 0
//<RM>: Added functions
// NOT being used //
template <class Image1T, class Image2T>
bool homography_ip_matching1( vw::ImageViewBase<Image1T> const& image1,
			       			  vw::ImageViewBase<Image2T> const& image2,
			       			  int ip_per_tile,
			       			  std::string const& output_name,
			       			  int inlier_threshold=10,
			       			  double nodata1 = std::numeric_limits<double>::quiet_NaN(),
			       			  double nodata2 = std::numeric_limits<double>::quiet_NaN(),
			       			  std::vector<ip::InterestPoint>& final_ip1 = NULL,
 			       			  std::vector<ip::InterestPoint>& final_ip2 = NULL );
Vector2i homography_rectification1( bool adjust_left_image_size,
			    					Vector2i const& left_size,
			    					Vector2i const& right_size,
			    					std::vector<ip::InterestPoint> const& left_ip,
			    					std::vector<ip::InterestPoint> const& right_ip,
			    					vw::Matrix<double>& left_matrix,
			    					vw::Matrix<double>& right_matrix,
			    					double threshRANSAC, 
									double minAvgDeltaY, 
									BBox2i bbox );
double calcAverageDeltaY( std::vector<Vector3> const& left_points, 
						  std::vector<Vector3> const& right_points);
Vector2i affine_epipolar_rectification1( Vector2i const& left_size,
                                 		 Vector2i const& right_size,
                                 		 std::vector<ip::InterestPoint> const& ip1,
                                 		 std::vector<ip::InterestPoint> const& ip2,
                                 		 Matrix<double>& left_matrix,
                                 		 Matrix<double>& right_matrix );
bool check_homography_matrix( Matrix<double>       const& left_matrix,
							  Matrix<double>       const& right_matrix,
			       			  std::vector<Vector3> const& left_points,
			       			  std::vector<Vector3> const& right_points,
			       			  double minAvgDeltaY, 
							  BBox2i bbox, double *avgDeltaY);
bool check_homography_matrix(Matrix<double>       const& H,
                               std::vector<Vector3> const& left_points,
                               std::vector<Vector3> const& right_points,
                               std::vector<size_t>  const& indices,
                               double minAvgDeltaY, BBox2i bbox);
double calcAverageDeltaY( std::vector<ip::InterestPoint> const& left_points, 
						  std::vector<ip::InterestPoint> const& right_points);
BBox2f calcSearchRange(std::vector<ip::InterestPoint> const& left_ip,
					   std::vector<ip::InterestPoint> const& right_ip, 
					   Matrix<double> const& left_matrix, 
					   Matrix<double> const& right_matrix, 
					   double multi);
BBox2f piecewiseAlignment_homography(ImageView<float> tile_left_image,
									 ImageView<float> tile_right_image,	
									 BBox2i bbox,
									 Vector2i& left_size,
									 Vector2i& right_size,
									 vw::Matrix<double>& left_matrix,
									 vw::Matrix<double>& right_matrix,
									 BBox2f local_search_range);
BBox2f piecewiseAlignment_affineepipolar( ImageView<float> tile_left_image,
										  ImageView<float> tile_right_image,	
										  BBox2i bbox,
										  Vector2i& left_size,
										  Vector2i& right_size,
										  vw::Matrix<double>& left_matrix,
										  vw::Matrix<double>& right_matrix,
										  BBox2f local_search_range);

/// Returns the properly cast cost mode type
stereo::CostFunctionType get_cost_mode_value() {
  switch(stereo_settings().cost_mode) {
    case 0: return stereo::ABSOLUTE_DIFFERENCE;
    case 1: return stereo::SQUARED_DIFFERENCE;
    case 2: return stereo::CROSS_CORRELATION;
    case 3: return stereo::CENSUS_TRANSFORM;
    case 4: return stereo::TERNARY_CENSUS_TRANSFORM;
    default: 
      vw_throw( ArgumentErr() << "Unknown value " << stereo_settings().cost_mode << " for cost-mode.\n" );
  };
}

/// Determine the proper subpixel mode to be used with SGM correlation
SemiGlobalMatcher::SgmSubpixelMode get_sgm_subpixel_mode() {

  switch(stereo_settings().subpixel_mode) {
    case  6: return SemiGlobalMatcher::SUBPIXEL_LINEAR;
    case  7: return SemiGlobalMatcher::SUBPIXEL_POLY4;
    case  8: return SemiGlobalMatcher::SUBPIXEL_COSINE;
    case  9: return SemiGlobalMatcher::SUBPIXEL_PARABOLA;
    case 10: return SemiGlobalMatcher::SUBPIXEL_NONE;
    case 11: return SemiGlobalMatcher::SUBPIXEL_LC_BLEND;
    default: return SemiGlobalMatcher::SUBPIXEL_LC_BLEND;
  };
}


// Read the search range from D_sub, and scale it to the full image
void read_search_range_from_dsub(ASPGlobalOptions & opt){

  // No D_sub is generated or should be used for seed mode 0.
  if (stereo_settings().seed_mode == 0)
    return;

  DiskImageView<vw::uint8> Lmask(opt.out_prefix + "-lMask.tif"),
                           Rmask(opt.out_prefix + "-rMask.tif");

  DiskImageView<PixelGray<float> > left_sub ( opt.out_prefix+"-L_sub.tif" ),
                                   right_sub( opt.out_prefix+"-R_sub.tif" );

  Vector2 downsample_scale( double(left_sub.cols()) / double(Lmask.cols()),
                            double(left_sub.rows()) / double(Lmask.rows()) );

  std::string d_sub_file = opt.out_prefix + "-D_sub.tif";
  if (!fs::exists(d_sub_file))
    return;

  ImageView<PixelMask<Vector2f> > sub_disp;
  read_image(sub_disp, d_sub_file);
  BBox2i search_range = stereo::get_disparity_range( sub_disp );
  search_range.min() = floor(elem_quot(search_range.min(),downsample_scale));
  search_range.max() = ceil (elem_quot(search_range.max(),downsample_scale));
  stereo_settings().search_range = search_range;
  
  vw_out() << "\t--> Read search range from D_sub: " << search_range << "\n";
}



/// Produces the low-resolution disparity file D_sub
void produce_lowres_disparity( ASPGlobalOptions & opt ) {

  // Set up handles to read the input images
  DiskImageView<vw::uint8> Lmask(opt.out_prefix + "-lMask.tif"),
                           Rmask(opt.out_prefix + "-rMask.tif");

  DiskImageView<PixelGray<float> > left_sub ( opt.out_prefix+"-L_sub.tif" ),
                                   right_sub( opt.out_prefix+"-R_sub.tif" );

  DiskImageView<uint8> left_mask_sub ( opt.out_prefix+"-lMask_sub.tif" ),
                       right_mask_sub( opt.out_prefix+"-rMask_sub.tif" );

  Vector2 downsample_scale( double(left_sub.cols()) / double(Lmask.cols()),
                            double(left_sub.rows()) / double(Lmask.rows()) );
  double mean_scale = (downsample_scale[0] + downsample_scale[1]) / 2.0;

  // Compute the initial search range in the subsampled image
  BBox2i search_range( floor(elem_prod(downsample_scale,stereo_settings().search_range.min())),
                       ceil (elem_prod(downsample_scale,stereo_settings().search_range.max())) );

  if ( stereo_settings().seed_mode == 1 ) {

    // Use low-res correlation to get the low-res disparity
    Vector2i expansion( search_range.width(),
                  			search_range.height() );
    expansion *= stereo_settings().seed_percent_pad / 2.0f;
    // Expand by the user selected amount. Default is 25%.
    search_range.min() -= expansion;
    search_range.max() += expansion;
    //VW_OUT(DebugMessage,"asp") << "D_sub search range: " << search_range << " px\n";
    vw_out() << "D_sub search range: " << search_range << " px\n";
    stereo::CostFunctionType cost_mode = get_cost_mode_value();
    Vector2i kernel_size  = stereo_settings().corr_kernel;
    int corr_timeout      = 5*stereo_settings().corr_timeout; // 5x, so try hard
    const int rm_half_kernel = 5; // Filter kernel size used by CorrelationView
    double seconds_per_op = 0.0;
    if (corr_timeout > 0)
      seconds_per_op = calc_seconds_per_op(cost_mode, left_sub, right_sub, kernel_size);

    SemiGlobalMatcher::SgmSubpixelMode sgm_subpixel_mode = get_sgm_subpixel_mode();
    Vector2i sgm_search_buffer = stereo_settings().sgm_search_buffer;;

    if (stereo_settings().rm_quantile_multiple <= 0.0)
    {
      // If we can process the entire image in one tile, don't use a collar.
      int collar_size = stereo_settings().sgm_collar_size;
      if ((opt.raster_tile_size[0] > left_sub.cols()) &&
          (opt.raster_tile_size[1] > left_sub.rows())   )
        collar_size = 0;
    
      // Warning: A giant function call approaches!
      // TODO: Why the extra filtering step here? PyramidCorrelationView already performs 1-3 iterations of outlier removal.
      std::string d_sub_file = opt.out_prefix + "-D_sub.tif";
      vw_out() << "Writing: " << d_sub_file << std::endl;
      vw::cartography::block_write_gdal_image( // Write to disk
          d_sub_file,
          rm_outliers_using_thresh( // Throw out individual pixels that are far from any neighbors
              vw::stereo::pyramid_correlate( // Compute image correlation using the PyramidCorrelationView class
                  left_sub, right_sub,
                  left_mask_sub, right_mask_sub,
                  vw::stereo::PREFILTER_LOG, stereo_settings().slogW,
                  search_range, kernel_size, cost_mode,
                  corr_timeout, seconds_per_op,
                  stereo_settings().xcorr_threshold, 
                  stereo_settings().min_xcorr_level,
                  rm_half_kernel,
                  stereo_settings().corr_max_levels,
                  static_cast<vw::stereo::CorrelationAlgorithm>(stereo_settings().stereo_algorithm),
                  collar_size, sgm_subpixel_mode, sgm_search_buffer, stereo_settings().corr_memory_limit_mb,
                  stereo_settings().corr_blob_filter_area*mean_scale,
                  stereo_settings().stereo_debug
              ),
              // To do: all these hard-coded values must be replaced with
              // appropriate params from user's stereo.default, for
              // consistency with how disparity is filtered in stereo_fltr,
              // when invoking disparity_cleanup_using_thresh.
              1, 1, // in stereo.default we have 5 5
              // Changing below the hard-coded value from 2.0 to using a
              // param.  The default value will still be 2.0 but is now
              // modifiable. Need to get rid of the 2.0/3.0 factor and
              // study how it affects the result.
              stereo_settings().rm_threshold*2.0/3.0,
              // Another change of hard-coded value to param. Get rid of 0.5/0.6
              // and study the effect.
              (stereo_settings().rm_min_matches/100.0)*0.5/0.6
          ), // End outlier removal arguments
          opt,
          TerminalProgressCallback("asp", "\t--> Low-resolution disparity:")
      );
      // End of giant function call block
    }
    else { // Use quantile based filtering - This filter needs to be profiled to improve its speed.
    
      // Compute image correlation using the PyramidCorrelationView class
      ImageView< PixelMask<Vector2f> > disp_image = vw::stereo::pyramid_correlate( 
                  left_sub, right_sub,
                  left_mask_sub, right_mask_sub,
                  vw::stereo::PREFILTER_LOG, stereo_settings().slogW,
                  search_range, kernel_size, cost_mode,
                  corr_timeout, seconds_per_op,
                  stereo_settings().xcorr_threshold, 
                  stereo_settings().min_xcorr_level,
                  rm_half_kernel,
                  stereo_settings().corr_max_levels,
                  static_cast<vw::stereo::CorrelationAlgorithm>(stereo_settings().stereo_algorithm), 
                  0, // No collar here, the entire image is written at once.
                  sgm_subpixel_mode, sgm_search_buffer, stereo_settings().corr_memory_limit_mb,
                  0, // Don't combine blob filtering with quantile filtering
                  stereo_settings().stereo_debug
              );

      std::string d_sub_file = opt.out_prefix + "-D_sub.tif";
      vw_out() << "Writing: " << d_sub_file << std::endl;
      vw::cartography::write_gdal_image( // Write to disk while removing outliers
          d_sub_file,
          rm_outliers_using_quantiles( // Throw out individual pixels that are far from any neighbors
              disp_image,
              stereo_settings().rm_quantile_percentile, stereo_settings().rm_quantile_multiple
          ),
          opt,
          TerminalProgressCallback("asp", "\t--> Low-resolution disparity:")
      );
    }

  }else if ( stereo_settings().seed_mode == 2 ) {
    // Use a DEM to get the low-res disparity
    boost::shared_ptr<camera::CameraModel> left_camera_model, right_camera_model;
    opt.session->camera_models(left_camera_model, right_camera_model);
    produce_dem_disparity(opt, left_camera_model, right_camera_model, opt.session->name());
  }else if ( stereo_settings().seed_mode == 3 ) {
    // D_sub is already generated by now by sparse_disp
  }

  read_search_range_from_dsub(opt); // TODO: We already call this when needed!
} // End produce_lowres_disparity


/// Adjust IP lists if alignment matrices are present.
double adjust_ip_for_align_matrix(std::string               const& out_prefix,
                                  vector<ip::InterestPoint>      & ip_left,
                                  vector<ip::InterestPoint>      & ip_right,
                                  double                    const  ip_scale) {

  // Check for alignment files
  bool left_align  = fs::exists(out_prefix+"-align-L.exr");
  bool right_align = fs::exists(out_prefix+"-align-R.exr");
  if (!left_align && !right_align)
    return ip_scale; // No alignment files -> Nothing to do.

  // Load alignment matrices
  Matrix<double> align_left_matrix  = math::identity_matrix<3>();
  Matrix<double> align_right_matrix = math::identity_matrix<3>();
  if (left_align)
    read_matrix(align_left_matrix, out_prefix + "-align-L.exr");
  if (right_align)
    read_matrix(align_right_matrix, out_prefix + "-align-R.exr");

  // Loop through all the IP we found
  for ( size_t i = 0; i < ip_left.size(); i++ ) {
    // Apply the alignment transforms to the recorded IP
    Vector3 l = align_left_matrix  * Vector3(ip_left [i].x, ip_left [i].y, 1);
    Vector3 r = align_right_matrix * Vector3(ip_right[i].x, ip_right[i].y, 1);

    // Normalize the coordinates, but don't divide by 0
    if (l[2] == 0 || r[2] == 0) 
      continue;
    l /= l[2];
    r /= r[2];

    ip_left [i].x = l[0];
    ip_left [i].y = l[1];
    ip_left [i].ix = l[0];
    ip_left [i].iy = l[1];
    
    ip_right[i].x = r[0];
    ip_right[i].y = r[1];
    ip_right[i].ix = r[0];
    ip_right[i].iy = r[1];
  }
  return 1.0; // If alignment files are present they take care of the scaling.
	      
} // End adjust_ip_for_align_matrix


/// Adjust IP lists if epipolar alignment was applied after the IP were created.
/// - Currently this condition can only happen if an IP file is inserted into the run
///   folder from another source such as bundle adjust!
/// - Returns true if any change was made to the interest points.
bool adjust_ip_for_epipolar_transform(ASPGlobalOptions          const& opt,
                                      std::string               const& match_file,
                                      vector<ip::InterestPoint>      & ip_left,
                                      vector<ip::InterestPoint>      & ip_right) {

  bool usePinholeEpipolar = ( (stereo_settings().alignment_method == "epipolar") &&
                              ( opt.session->name() == "pinhole" ||
                                opt.session->name() == "nadirpinhole") );
  
  if (!usePinholeEpipolar) 
    return false;
  
  // This function does nothing if we are not using epipolar alignment,
  //  or if the IP were found using one of the aligned images.
  const std::string sub_match_file     = opt.out_prefix + "-L_sub__R_sub.match";
  const std::string aligned_match_file = opt.out_prefix + "-L__R.match";
  if ( (stereo_settings().alignment_method != "epipolar") ||
       (match_file == sub_match_file) || (match_file == aligned_match_file) )
    return false;

  vw_out() << "Applying epipolar adjustment to input IP match file...\n";

  // Need to cast the session pointer to Pinhole type to access the function we need.
  StereoSessionPinhole* pinPtr = dynamic_cast<StereoSessionPinhole*>(opt.session.get());
  if (pinPtr == NULL) 
    vw_throw(ArgumentErr() << "Expected a pinhole camera.\n");
  // Must initialize below the two cameras to something to respect the constructor.
  asp::PinholeCamTrans trans_left = asp::PinholeCamTrans(vw::camera::PinholeModel(),
                                                         vw::camera::PinholeModel());
  asp::PinholeCamTrans trans_right = trans_left;
  pinPtr->pinhole_cam_trans(trans_left, trans_right);

  // Apply the transforms to all the IP we found
  for ( size_t i = 0; i < ip_left.size(); i++ ) {

    Vector2 ip_in_left (ip_left [i].x, ip_left [i].y);
    Vector2 ip_in_right(ip_right[i].x, ip_right[i].y);

    Vector2 ip_out_left  = trans_left.forward(ip_in_left);
    Vector2 ip_out_right = trans_right.forward(ip_in_right);

    ip_left [i].x = ip_out_left [0]; // Store transformed points
    ip_left [i].y = ip_out_left [1];
    ip_right[i].x = ip_out_right[0];
    ip_right[i].y = ip_out_right[1];
  }

  return true;
} // End adjust_ip_for_epipolar_transform



// TODO: Duplicate of hidden function in vw/src/InterestPoint/Matcher.cc!
std::string strip_path(std::string out_prefix, std::string filename){

  // If filename starts with out_prefix followed by dash, strip both.
  // Also strip filename extension.

  std::string ss = out_prefix + "-";
  size_t found = filename.find(ss);

  if (found != std::string::npos)
    filename.erase(found, ss.length());

  filename = fs::path(filename).stem().string();

  return filename;
}

/// Detect IP in the _sub images or the original images if they are not too large.
/// - Usually an IP file is written in stereo_pprc, but for some input scenarios
///   this function will need to be used to generate them here.
/// - The input match file path can be changed depending on what exists on disk.
/// - Returns the scale from the image used for IP to the full size image.
/// - The binary interest point file will be written to disk.
double compute_ip(ASPGlobalOptions & opt, std::string & match_filename) {

  vw_out() << "\t    * Loading images for IP detection.\n";

  // Choose whether to use the full or _sub images

  // Use the full image if all dimensions are smaller than this.
  const int SIZE_CUTOFF = 8000;

  const std::string left_image_path_full  = opt.out_prefix+"-L.tif";
  const std::string right_image_path_full = opt.out_prefix+"-R.tif";
  const std::string left_image_path_sub   = opt.out_prefix+"-L_sub.tif";
  const std::string right_image_path_sub  = opt.out_prefix+"-R_sub.tif";

  // TODO: Just call the right function everywhere rather than computing its result by hand.
  const std::string full_match_file       = ip::match_filename(opt.out_prefix, opt.in_file1, opt.in_file2);
  const std::string sub_match_file        = opt.out_prefix + "-L_sub__R_sub.match";
  const std::string aligned_match_file    = opt.out_prefix + "-L__R.match";

  // TODO: The logic below is wrong. Don't read the first match file
  // that happens to exist on disk and hope for the best.  That could
  // be an incorrect one. At this stage we know exactly the files that
  // need processing. Check if the desired file exists, and read that
  // one, or create it if missing.
  
  // Try the full match file first
  if (fs::exists(full_match_file)) {
    vw_out() << "IP file found: " << full_match_file << std::endl;
    match_filename = full_match_file;
    return 1.0;
  }

  // TODO: Unify with function in vw/src/InterestPoint/Matcher.h!
  // filenames longer than this must be chopped, as too long names
  // cause problems later with boost.
  int max_len = 40;
  std::string name1 = strip_path(opt.out_prefix, opt.in_file1).substr(0, max_len);
  std::string name2 = strip_path(opt.out_prefix, opt.in_file2).substr(0, max_len);

  // Next try the cropped match file names which will be at full scale.
  // TODO: This is unnecessary. Just call the right function to find
  // the match file.
  std::vector<std::string> match_names;
  match_names.push_back(opt.out_prefix + "-L-cropped__R-cropped.match");
  match_names.push_back(opt.out_prefix + "-"+name1+"__R-cropped.match");
  match_names.push_back(opt.out_prefix + "-L-cropped__"+name2+".match");
  match_names.push_back(aligned_match_file);
  for (size_t i=0; i<match_names.size(); ++i) {
    if (fs::exists(match_names[i])) {
      vw_out() << "IP file found: " << match_names[i] << std::endl;
      match_filename = match_names[i];
      return 1.0;
    }
  }

  // Now try the sub match file, which requires us to compute the scale.
  std::string left_image_path  = left_image_path_full;
  std::string right_image_path = right_image_path_full;
  Vector2i full_size     = file_image_size(left_image_path_full);
  bool     use_full_size = (((full_size[0] < SIZE_CUTOFF) && (full_size[1] < SIZE_CUTOFF))
                            || ((stereo_settings().alignment_method != "epipolar") &&
                                (stereo_settings().alignment_method != "none"    )   ));
  // Other alignment methods find IP in the stereo_pprc phase using the full size.

  // Compute the scale.
  double ip_scale = 1.0;
  if (!use_full_size) {
    left_image_path  = left_image_path_sub;
    right_image_path = right_image_path_sub;
    
    ip_scale = sum(elem_quot( Vector2(file_image_size( opt.out_prefix+"-L_sub.tif" )),
                              Vector2(file_image_size( opt.out_prefix+"-L.tif" ) ) )) +
               sum(elem_quot( Vector2(file_image_size( opt.out_prefix+"-R_sub.tif" )),
                              Vector2(file_image_size( opt.out_prefix+"-R.tif" ) ) ));
    ip_scale /= 4.0f;
    match_filename = sub_match_file; // If not using full size we should expect this file
    
    // Check for the file.
    if (fs::exists(sub_match_file)) {
      vw_out() << "IP file found: " << sub_match_file << std::endl;
      return ip_scale;
    }
  }
  else
    match_filename = aligned_match_file;

  vw_out() << "No IP file found, computing IP now.\n";
  
  // Load the images
  boost::shared_ptr<DiskImageResource> left_rsrc (DiskImageResourcePtr(left_image_path )),
                                       right_rsrc(DiskImageResourcePtr(right_image_path));

  // Read the no-data values written to disk previously when
  // the normalized left and right sub-images were created.
  float left_nodata_value  = numeric_limits<float>::quiet_NaN();
  float right_nodata_value = numeric_limits<float>::quiet_NaN();
  if (left_rsrc->has_nodata_read ()) left_nodata_value  = left_rsrc->nodata_read();
  if (right_rsrc->has_nodata_read()) right_nodata_value = right_rsrc->nodata_read();
  
  // These images should be small enough to fit in memory
  ImageView<float> left_image  = DiskImageView<float>(left_rsrc);
  ImageView<float> right_image = DiskImageView<float>(right_rsrc);

  // No interest point operations have been performed before
  vw_out() << "\t    * Locating Interest Points\n";

  // Use this code in a relatively specific case
  // - Only tested with IceBridge data so far!
  // - Some changes will be required for this to work in more general cases.
  bool success = false;
  if (use_full_size && opt.session->is_nadir_facing() && 
      (stereo_settings().alignment_method == "epipolar") ) {

    // Load camera models
    boost::shared_ptr<camera::CameraModel> left_camera_model, right_camera_model;
    opt.session->camera_models(left_camera_model, right_camera_model);
    
    // Obtain the datum
    const bool use_sphere_for_isis = false;
    cartography::Datum datum = opt.session->get_datum(left_camera_model.get(), use_sphere_for_isis);

    // Since these are epipolar aligned images it should be small
    double epipolar_threshold = 5;
    if (stereo_settings().epipolar_threshold > 0)
      epipolar_threshold = stereo_settings().epipolar_threshold;

    const bool single_threaded_camera = false;
    success = ip_matching_no_align(single_threaded_camera,
                                   left_camera_model.get(), right_camera_model.get(),
                                   left_image, right_image,
                                   stereo_settings().ip_per_tile,
                                   datum, match_filename, epipolar_threshold,
                                   stereo_settings().ip_uniqueness_thresh,
                                   left_nodata_value, right_nodata_value);
  } // End nadir epipolar full image case
  else {
    // In all other cases, run a more general IP matcher.
    
    // TODO: Depending on alignment method, we can tailor the IP filtering strategy.
    double thresh_factor = stereo_settings().ip_inlier_factor; // 1/15 by default

    // This range is extra large to handle elevation differences.
    const int inlier_threshold = 200*(15.0*thresh_factor);  // 200 by default
    
    success = asp::homography_ip_matching(left_image, right_image,
                                          stereo_settings().ip_per_tile,
                                          match_filename, inlier_threshold,
                                          left_nodata_value, right_nodata_value);
  }

  if (!success)
    vw_throw(ArgumentErr() << "Could not find interest points.\n");

  return ip_scale;
}



BBox2i get_search_range_from_ip_hists(vw::math::Histogram const& hist_x,
                                      vw::math::Histogram const& hist_y,
                                      double edge_discard_percentile = 0.05) {

  const double min_percentile = edge_discard_percentile;
  const double max_percentile = 1.0 - edge_discard_percentile;

  const Vector2 FORCED_EXPANSION = Vector2(30,2); // Must expand range by at least this much
  double search_scale = 2.0;
  size_t min_bin_x = hist_x.get_percentile(min_percentile);
  size_t min_bin_y = hist_y.get_percentile(min_percentile);
  size_t max_bin_x = hist_x.get_percentile(max_percentile);
  size_t max_bin_y = hist_y.get_percentile(max_percentile);
  Vector2 search_min(hist_x.get_bin_center(min_bin_x), hist_y.get_bin_center(min_bin_y));
  Vector2 search_max(hist_x.get_bin_center(max_bin_x), hist_y.get_bin_center(max_bin_y));
  Vector2 search_center = (search_max + search_min) / 2.0;
  Vector2 d_min = search_min - search_center; // TODO: Make into a bbox function!
  Vector2 d_max = search_max - search_center;
  
  vw_out(InfoMessage,"asp") << "Percentile filtered range: " 
                            << BBox2i(d_min, d_max) << std::endl;
  // Enforce a minimum expansion on the search range in each direction
  Vector2 min_expand = d_min*search_scale;
  Vector2 max_expand = d_max*search_scale;
  for (int i=0; i<2; ++i) {
    if (min_expand[i] > -1*FORCED_EXPANSION[i])
      min_expand[i] = -1*FORCED_EXPANSION[i];
    if (max_expand[i] < FORCED_EXPANSION[i])
      max_expand[i] = FORCED_EXPANSION[i];
  }
  
  search_min = search_center + min_expand;
  search_max = search_center + max_expand;
  Vector2i search_minI(floor(search_min[0]), floor(search_min[1])); // Round outwards
  Vector2i search_maxI(ceil (search_max[0]), ceil (search_max[1]));
/*
   // Debug code to print all the points
  for (size_t i = 0; i < matched_ip1.size(); i++) {
    Vector2f diff(i_scale * (matched_ip2[i].x - matched_ip1[i].x), 
                  i_scale * (matched_ip2[i].y - matched_ip1[i].y));
    //Vector2f diff(matched_ip2[i].x - matched_ip1[i].x, 
    //              matched_ip2[i].y - matched_ip1[i].y);
    vw_out(InfoMessage,"asp") << matched_ip1[i].x <<", "<<matched_ip1[i].y 
              << " <> " 
              << matched_ip2[i].x <<", "<<matched_ip2[i].y 
               << " DIFF " << diff << endl;
  }
*/
  
  //vw_out(InfoMessage,"asp") << "i_scale is : "       << i_scale << endl;
  
  return BBox2i(search_minI, search_maxI);
}
  


/// Use existing interest points to compute a search range
/// - This function could use improvement!
/// - Should it be used in all cases?
BBox2i approximate_search_range(ASPGlobalOptions & opt, 
                                double ip_scale, std::string const& match_filename) {

  vw_out() << "\t--> Using interest points to determine search window.\n";
  vector<ip::InterestPoint> in_ip1, in_ip2, matched_ip1, matched_ip2;

  // The interest points must have been created outside this function
  if (!fs::exists(match_filename))
    vw_throw( ArgumentErr() << "Missing IP file: " << match_filename);

  vw_out() << "\t    * Loading match file: " << match_filename << "\n";
  ip::read_binary_match_file(match_filename, in_ip1, in_ip2);

  // TODO: Consolidate IP adjustment
  // TODO: This logic is messed up. We __know__ from stereo_settings() what
  // alignment method is being used and what scale we are at, there is no
  // need to try to read various and likely old files from disk
  // to infer that. You can get the wrong answer.
  
  // Handle alignment matrices if they are present
  // - Scale is reset to 1.0 if alignment matrices are present.
  ip_scale = adjust_ip_for_align_matrix(opt.out_prefix, in_ip1, in_ip2, ip_scale);
  vw_out() << "\t    * IP computed at scale: " << ip_scale << "\n";
  float i_scale = 1.0/ip_scale;

  // Adjust the IP if they came from input images and these images are epipolar aligned
  adjust_ip_for_epipolar_transform(opt, match_filename, in_ip1, in_ip2);

  // Filter out IPs which fall outside the specified elevation range
  boost::shared_ptr<camera::CameraModel> left_camera_model, right_camera_model;
  opt.session->camera_models(left_camera_model, right_camera_model);
  cartography::Datum datum = opt.session->get_datum(left_camera_model.get(), false);

  // We already corrected for align matrix, so transforms should be identity here.
  vw::TransformRef left_tx  = vw::TransformRef(vw::TranslateTransform(0,0));
  vw::TransformRef right_tx = vw::TransformRef(vw::TranslateTransform(0,0));

  // Filter out IPs which fall outside the specified elevation and lonlat range
  // TODO: Don't do this with cropped input images!!!!!
  size_t num_left = asp::filter_ip_by_lonlat_and_elevation(left_camera_model.get(),
							   right_camera_model.get(),
							   datum, in_ip1, in_ip2,
							   left_tx, right_tx, ip_scale,
							   stereo_settings().elevation_limit,
							   stereo_settings().lon_lat_limit,
							   matched_ip1, matched_ip2);

  // If the user set this, filter by disparity of ip.
  // TODO: This kind of logic is present below one more time, at
  // get_search_range_from_ip_hists() where a factor of 2 is used!
  // This logic better be integrated together!
  Vector2 disp_params = stereo_settings().remove_outliers_by_disp_params;
  if (disp_params[0] < 100.0) // not enabled by default
    asp::filter_ip_by_disparity(disp_params[0], disp_params[1], matched_ip1, matched_ip2); 
  
  // Quit if we don't have the requested number of IP.
  if (static_cast<int>(num_left) < stereo_settings().min_num_ip)
    vw_throw(ArgumentErr() << "Number of IPs left after filtering is " << num_left
                           << " which is less than the required amount of " 
                           << stereo_settings().min_num_ip << ", aborting stereo_corr.\n");

  // Find search window based on interest point matches
  size_t num_ip = matched_ip1.size();
  vw_out(InfoMessage,"asp") << "Estimating search range with: " 
                            << num_ip << " interest points.\n";

  // Record the disparities for each point pair
  const double BIG_NUM   =  99999999;
  const double SMALL_NUM = -99999999;
  std::vector<double> dx, dy;
  double min_dx = BIG_NUM, max_dx = SMALL_NUM,
         min_dy = BIG_NUM, max_dy = SMALL_NUM;
  for (size_t i = 0; i < num_ip; i++) {
    double diffX = i_scale * (matched_ip2[i].x - matched_ip1[i].x);
    double diffY = i_scale * (matched_ip2[i].y - matched_ip1[i].y);      
    dx.push_back(diffX);
    dy.push_back(diffY);
    if (diffX < min_dx) min_dx = diffX;
    if (diffY < min_dy) min_dy = diffY;
    if (diffX > max_dx) max_dx = diffX;
    if (diffY > max_dy) max_dy = diffY;
  }

  vw_out(InfoMessage,"asp") << "Initial search range: " 
        << BBox2i(Vector2(min_dx,min_dy),Vector2(max_dx,max_dy)) << std::endl;

  const int MAX_SEARCH_WIDTH = 4000; // Try to avoid searching this width
  const int MIN_SEARCH_WIDTH = 200;  // Under this width don't filter IP.
  const Vector2i MINIMAL_EXPAND(10,1);

  // If the input search range is small just expand it a bit and
  //  return without doing any filtering.  
  if (max_dx-min_dx <= MIN_SEARCH_WIDTH) {
    BBox2i search_range(Vector2i(min_dx,min_dy),Vector2i(max_dx,max_dy));
    search_range.min() -= MINIMAL_EXPAND; // BBox2.expand() function does not always work!!!!
    search_range.max() += MINIMAL_EXPAND;
    vw_out(InfoMessage,"asp") << "Using expanded search range: " 
                              << search_range << std::endl;
    return search_range;
  }
  
  // Compute histograms
  const int NUM_BINS = 2000; // Accuracy is important with scaled pixels
  vw::math::Histogram hist_x(NUM_BINS, min_dx, max_dx);
  vw::math::Histogram hist_y(NUM_BINS, min_dy, max_dy);
  for (size_t i = 0; i < dx.size(); ++i){
    hist_x.add_value(dx[i]);
    hist_y.add_value(dy[i]);
  }


  //printf("min x,y = %lf, %lf, max x,y = %lf, %lf\n", min_dx, min_dy, max_dx, max_dy);
  
  //for (int i=0; i<NUM_BINS; ++i) {
  //  printf("%d => X: %lf: %lf,   Y:  %lf: %lf\n", i, centers_x[i], hist_x[i], centers_y[i], hist_y[i]);
  //}

  const double PERCENTILE_CUTOFF     = 0.05; // Gradually increase the filtering
  const double PERCENTILE_CUTOFF_INC = 0.05; //  until the search width is reasonable.
  const double MAX_PERCENTILE_CUTOFF = 0.201;

  double current_percentile_cutoff = PERCENTILE_CUTOFF;
  int search_width = MAX_SEARCH_WIDTH + 1;
  BBox2i search_range;
  while (true) {
    vw_out() << "Filtering IP with percentile cutoff " << current_percentile_cutoff << std::endl;
    search_range = get_search_range_from_ip_hists(hist_x, hist_y, current_percentile_cutoff);
    vw_out() << "Scaled search range = " << search_range << std::endl;
    search_width = search_range.width();
    
    // Increase the percentile cutoff in case we need to filter out more IP
    current_percentile_cutoff += PERCENTILE_CUTOFF_INC;
    if (current_percentile_cutoff > MAX_PERCENTILE_CUTOFF) {
      if (search_width < MAX_SEARCH_WIDTH)
        vw_out() << "Exceeded maximum filter cutoff of " << MAX_PERCENTILE_CUTOFF
                 << ", keeping current search range\n";
      break; // No more filtering is possible, exit the loop.
    }
      
    if (search_width < MAX_SEARCH_WIDTH)
      break; // Happy with search range, exit the loop.
    else
      vw_out() << "Search width of " << search_width << " is greater than desired limit of "
               << MAX_SEARCH_WIDTH << ", retrying with more aggressive IP filter\n";
  } // End search range determination loop

  
  // Prevent any dimension from being length zero,
  //  otherwise future parts to ASP will fail.
  // TODO: Fix ASP and SGM handling of small boxes!
  //       - Currently code has a minimum search height of 5!
  if (search_range.empty())
    vw_throw(ArgumentErr() << "Computed an empty search range!");
  
  return search_range;
} // End function approximate_search_range


/// The first step of correlation computation.
void lowres_correlation( ASPGlobalOptions & opt ) {

  vw_out() << "\n[ " << current_posix_time_string() << " ] : Stage 1 --> LOW-RESOLUTION CORRELATION \n";

  // Working out search range if need be
  if (stereo_settings().is_search_defined()) {
    vw_out() << "\t--> Using user-defined search range.\n";

    // Update user provided search range based on input crops
    bool crop_left  = (stereo_settings().left_image_crop_win  != BBox2i(0, 0, 0, 0));
    bool crop_right = (stereo_settings().right_image_crop_win != BBox2i(0, 0, 0, 0));
    if (crop_left && !crop_right)
      stereo_settings().search_range += stereo_settings().left_image_crop_win.min();
    if (!crop_left && crop_right)
      stereo_settings().search_range -= stereo_settings().right_image_crop_win.min();

  }else if (stereo_settings().seed_mode == 2){
    // Do nothing as we will compute the search range based on D_sub
  }else if (stereo_settings().seed_mode == 3){
    // Do nothing as low-res disparity (D_sub) is already provided by sparse_disp
  } else { // Regular seed mode

    // If there is no match file for the input images, gather some IP from the
    // low resolution images. This routine should only run for:
    //   Pinhole + Epipolar
    //   Alignment method none
    //   Cases where either input image is cropped, in which case the IP name is different.
    // Everything else should gather IP's all the time during stereo_pprc.
    // - TODO: When inputs are cropped, use the cropped IP!
    

    // Compute new IP and write them to disk.
    // - If IP are already on disk this function will load them instead.
    // - This function will choose an appropriate IP computation based on the input images.
    string match_filename;
    double ip_scale;
    ip_scale = compute_ip(opt, match_filename);

    // This function applies filtering to find good points
    stereo_settings().search_range = approximate_search_range(opt, ip_scale, match_filename);
  
 
    vw_out() << "\t--> Detected search range: " << stereo_settings().search_range << "\n";
  } // End of case where we had to calculate the search range

  // If the user specified a search range limit, apply it here.
  if ((stereo_settings().search_range_limit.min() != Vector2i()) || 
      (stereo_settings().search_range_limit.max() != Vector2i())   ) {     
    stereo_settings().search_range.crop(stereo_settings().search_range_limit);
    vw_out() << "\t--> Detected search range constrained to: " << stereo_settings().search_range << "\n";
  }

  // At this point stereo_settings().search_range is populated

  DiskImageView<vw::uint8> Lmask(opt.out_prefix + "-lMask.tif"),
                           Rmask(opt.out_prefix + "-rMask.tif");

  // Performing disparity on sub images
  if ( stereo_settings().seed_mode > 0 ) {

    // Reuse prior existing D_sub if it exists, unless we
    // are cropping the images each time, when D_sub must
    // be computed anew each time.
    bool crop_left  = (stereo_settings().left_image_crop_win  != BBox2i(0, 0, 0, 0));
    bool crop_right = (stereo_settings().right_image_crop_win != BBox2i(0, 0, 0, 0));
    bool rebuild    = crop_left || crop_right;

    string sub_disp_file = opt.out_prefix+"-D_sub.tif";
    try {
      vw_log().console_log().rule_set().add_rule(-1,"fileio");
      DiskImageView<PixelMask<Vector2f> > test(sub_disp_file);
      vw_settings().reload_config();
    } catch (vw::IOErr const& e) {
      vw_settings().reload_config();
      rebuild = true;
    } catch (vw::ArgumentErr const& e ) {
      // Throws on a corrupted file.
      vw_settings().reload_config();
      rebuild = true;
    }

    if ( rebuild )
      produce_lowres_disparity(opt); // Note: This does not always remake D_sub!
    else
      vw_out() << "\t--> Using cached low-resolution disparity: " << sub_disp_file << "\n";
  }

  vw_out() << "\n[ " << current_posix_time_string() << " ] : LOW-RESOLUTION CORRELATION FINISHED \n";
} // End lowres_correlation



/// This correlator takes a low resolution disparity image as an input
/// so that it may narrow its search range for each tile that is processed.

//<Ricardo Monteiro>: added m_local_hom_L and m_local_size as inputs to handle affineepipolar piecewise alignment
class SeededCorrelatorView : public ImageViewBase<SeededCorrelatorView> {
  DiskImageView<PixelGray<float> >   m_left_image;
  DiskImageView<PixelGray<float> >   m_right_image;
  DiskImageView<vw::uint8> m_left_mask;
  DiskImageView<vw::uint8> m_right_mask;
  ImageViewRef<PixelMask<Vector2f> > m_sub_disp;
  ImageViewRef<PixelMask<Vector2i> > m_sub_disp_spread;
  ImageView<Matrix3x3> & m_local_hom;
  ImageView<Matrix3x3> & m_local_hom_L;
  ImageView<Matrix3x3> & m_local_size;

  // Settings
  Vector2  m_upscale_factor;
  BBox2i   m_seed_bbox;
  Vector2i m_kernel_size;
  stereo::CostFunctionType m_cost_mode;
  int      m_corr_timeout;
  double   m_seconds_per_op;

public:

  // Set these input types here instead of making them template arguments
  typedef DiskImageView<PixelGray<float> >   ImageType;
  typedef DiskImageView<vw::uint8>           MaskType;
  typedef ImageViewRef<PixelMask<Vector2f> > DispSeedImageType;
  typedef ImageViewRef<PixelMask<Vector2i> > SpreadImageType;
  typedef ImageType::pixel_type InputPixelType;

  SeededCorrelatorView( ImageType             const& left_image,
                        ImageType             const& right_image,
                        MaskType              const& left_mask,
                        MaskType              const& right_mask,
                        DispSeedImageType     const& sub_disp,
                        SpreadImageType       const& sub_disp_spread,
                        ImageView<Matrix3x3>  & local_hom,
						ImageView<Matrix3x3>  & local_hom_L,
						ImageView<Matrix3x3>  & local_size,
                        Vector2i const& kernel_size,
                        stereo::CostFunctionType cost_mode,
                        int corr_timeout, double seconds_per_op):
    m_left_image(left_image.impl()), m_right_image(right_image.impl()),
    m_left_mask (left_mask.impl ()), m_right_mask (right_mask.impl ()),
    m_sub_disp(sub_disp.impl()), m_sub_disp_spread(sub_disp_spread.impl()),
    m_local_hom(local_hom),
	m_local_hom_L(local_hom_L), m_local_size(local_size),
    m_kernel_size(kernel_size),  m_cost_mode(cost_mode),
    m_corr_timeout(corr_timeout), m_seconds_per_op(seconds_per_op){ 
    m_upscale_factor[0] = double(m_left_image.cols()) / m_sub_disp.cols();
    m_upscale_factor[1] = double(m_left_image.rows()) / m_sub_disp.rows();
    m_seed_bbox = bounding_box( m_sub_disp );
  }

  // Image View interface
  typedef PixelMask<Vector2f> pixel_type;
  typedef pixel_type          result_type;
  typedef ProceduralPixelAccessor<SeededCorrelatorView> pixel_accessor;

  inline int32 cols  () const { return m_left_image.cols(); }
  inline int32 rows  () const { return m_left_image.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

  inline pixel_type operator()(double /*i*/, double /*j*/, int32 /*p*/ = 0) const {
    vw_throw(NoImplErr() << "SeededCorrelatorView::operator()(...) is not implemented");
    return pixel_type();
  }
  	/// Does the work
  	typedef CropView<ImageView<pixel_type> > prerasterize_type;
  	inline prerasterize_type prerasterize(BBox2i const& bbox) const {

    bool use_local_homography = stereo_settings().use_local_homography;
    Matrix<double> fullres_hom = math::identity_matrix<3>();
    ImageView<InputPixelType> right_trans_img;
    ImageView<vw::uint8     > right_trans_mask;
	//<RM>: Added variables
	//<RM>: increase bbox where the piecewise alignment is applied to reduce border artifacts
	BBox2f new_local_search_range;
	int ts = ASPGlobalOptions::corr_tile_size();
	int margin = ts * 0.4; //<RM>: This can be lower to improve speed
#if DEBUG_RM
  	cout << "start of tile " << bbox << endl;
	char outputName[30];
	cartography::GdalWriteOptions geo_opt;
	int W = bbox.min().x()/ts, H = bbox.min().y()/ts;
#endif
	//<RM>: the left image tile might be transformed as well depending on the alignment technique
    ImageView<InputPixelType> left_trans_img; 
    ImageView<vw::uint8     > left_trans_mask;
	//<RM>: new bbox that includes a margin of pixels 
	BBox2i newBBox = bbox;
	newBBox.expand(margin);
	newBBox.crop(bounding_box(m_left_image));
#if DEBUG_RM
	cout << "[tile(" << H << "," << W << " newBBox = " << newBBox << endl;
	cout << "[tile(" << H << "," << W << " bbox = " << bbox << endl;
#endif
	//<RM>: the alignment is going to be applied only to the tile itself 
	ImageView<PixelGray<float> > tile_right_image = crop(m_right_image.impl(), newBBox);
	ImageView<PixelGray<float> > tile_left_image = crop(m_left_image.impl(), newBBox);
	ImageView<vw::uint8> tile_right_image_mask = crop(m_right_mask.impl(), newBBox);
	ImageView<vw::uint8> tile_left_image_mask = crop(m_left_mask.impl(), newBBox);
	//<RM>: the piecewise alignment is might require a transformation for each tile individually
	Matrix<double>  align_left_matrix  = math::identity_matrix<3>(),
                   	align_right_matrix = math::identity_matrix<3>();

    bool do_round = true; // round integer disparities after transform
    // User strategies
    BBox2f local_search_range;
    if ( stereo_settings().seed_mode > 0 ) {

      // The low-res version of bbox
      BBox2i seed_bbox( elem_quot(bbox.min(), m_upscale_factor),
			                  elem_quot(bbox.max(), m_upscale_factor) );
      seed_bbox.expand(1);
      seed_bbox.crop( m_seed_bbox );
      // Get the disparity range in d_sub corresponding to this tile.
      VW_OUT(DebugMessage, "stereo") << "Getting disparity range for : " << seed_bbox << "\n";
      DispSeedImageType disparity_in_box = crop( m_sub_disp, seed_bbox );

      if (!use_local_homography){
        local_search_range = stereo::get_disparity_range( disparity_in_box );
      }else{ // use local homography
        int ts = ASPGlobalOptions::corr_tile_size();
        lowres_hom = m_local_hom(bbox.min().x()/ts, bbox.min().y()/ts);
        local_search_range = stereo::get_disparity_range
          (transform_disparities(do_round, seed_bbox,
			     lowres_hom, disparity_in_box));
      }

      bool has_sub_disp_spread = ( m_sub_disp_spread.cols() != 0 &&
			                             m_sub_disp_spread.rows() != 0 );
      // Sanity check: If m_sub_disp_spread was provided, it better have the same size as sub_disp.
      if ( has_sub_disp_spread &&
           m_sub_disp_spread.cols() != m_sub_disp.cols() &&
           m_sub_disp_spread.rows() != m_sub_disp.rows() ){
        vw_throw( ArgumentErr() << "stereo_corr: D_sub and D_sub_spread must have equal sizes.\n");
      }

      if (has_sub_disp_spread){
        // Expand the disparity range by m_sub_disp_spread.
        SpreadImageType spread_in_box = crop( m_sub_disp_spread, seed_bbox );

        if (!use_local_homography){
          BBox2f spread = stereo::get_disparity_range( spread_in_box );
          local_search_range.min() -= spread.max();
          local_search_range.max() += spread.max();
        }else{
           //TODO: Check!
          DispSeedImageType upper_disp = transform_disparities(do_round, seed_bbox, lowres_hom,
                                                               disparity_in_box + spread_in_box);
          DispSeedImageType lower_disp = transform_disparities(do_round, seed_bbox, lowres_hom,
                                                               disparity_in_box - spread_in_box);
          BBox2f upper_range = stereo::get_disparity_range(upper_disp);
          BBox2f lower_range = stereo::get_disparity_range(lower_disp);

          local_search_range = upper_range;
          local_search_range.grow(lower_range);
        } //endif use_local_homography
      } //endif has_sub_disp_spread

      // TODO: Check!
      if (use_local_homography){
        Vector3 upscale(     m_upscale_factor[0],     m_upscale_factor[1], 1 );
        Vector3 dnscale( 1.0/m_upscale_factor[0], 1.0/m_upscale_factor[1], 1 );
        fullres_hom = diagonal_matrix(upscale)*lowres_hom*diagonal_matrix(dnscale);
		//<RM>: from this point on the fullres_hom is going to be overwritten by the tranformation calculated by a different piecewise alignment technique
		local_search_range = stereo::get_disparity_range( disparity_in_box );		
		Vector2i left_size = newBBox.size(), right_size = newBBox.size();
	    new_local_search_range = piecewiseAlignment_affineepipolar(tile_left_image, tile_right_image, bbox, left_size, right_size, align_left_matrix, align_right_matrix, local_search_range);
		//new_local_search_range = piecewiseAlignment_homography(tile_left_image, tile_right_image, bbox, left_size, right_size, align_left_matrix, align_right_matrix, local_search_range); //<RM>: Alternative method - not very good
		right_size = left_size;
		fullres_hom = align_right_matrix;
		//<RM>: write tranformation matrices for both left and right tiles to file
		m_local_hom(bbox.min().x()/ts, bbox.min().y()/ts) = fullres_hom; 
		m_local_hom_L(bbox.min().x()/ts, bbox.min().y()/ts) = align_left_matrix;
		//<RM>: write the aligned tile size to file to be accessed later by stereo_rfne
		m_local_size(bbox.min().x()/ts, bbox.min().y()/ts)(0,0) = left_size.x();
        m_local_size(bbox.min().x()/ts, bbox.min().y()/ts)(0,1) = left_size.y();
		//<RM>: tranform left tile
		ImageView< PixelMask<InputPixelType> > left_trans_masked_img = transform (copy_mask( tile_left_image.impl(),
																		  create_mask(tile_left_image_mask.impl()) ),
	               														  HomographyTransform(align_left_matrix),
																		  left_size.x(), left_size.y()); 
        left_trans_img  = apply_mask(left_trans_masked_img);
        left_trans_mask = channel_cast_rescale<uint8>(select_channel(left_trans_masked_img, 1));
		//<RM>: tranform right tile
        ImageView< PixelMask<InputPixelType> > right_trans_masked_img = transform (copy_mask(tile_right_image.impl(),
					 													   create_mask(tile_right_image_mask.impl()) ),
	               														   HomographyTransform(fullres_hom),
																		   right_size.x(), right_size.y()); 
        right_trans_img  = apply_mask(right_trans_masked_img);
        right_trans_mask = channel_cast_rescale<uint8>(select_channel(right_trans_masked_img, 1));
#if DEBUG_RM
		cout << "[tile(" << H << "," << W << " left_size after piecewise alignment = " << left_size << endl;
		cout << "[tile(" << H << "," << W << " right_size after piecewise alignment = " << right_size << endl;
		cout << "[tile(" << H << "," << W << " local_search_range = " << local_search_range << endl;
		cout << "[tile(" << H << "," << W << ") " << fullres_hom << "]" << endl;
		cout << "[tile(" << H << "," << W << ") " << align_left_matrix << "]" << endl;
		sprintf(outputName, "tile_R_%d_%d.tif", H, W);
		block_write_gdal_image(outputName, tile_right_image, geo_opt);
		sprintf(outputName, "tile_L_%d_%d.tif", H, W);
		block_write_gdal_image(outputName, tile_left_image, geo_opt);
		sprintf(outputName, "piecewiseHomography_R_%d_%d.tif", H, W);
		block_write_gdal_image(outputName, right_trans_img, geo_opt);
		sprintf(outputName, "piecewiseHomography_L_%d_%d.tif", H, W);
		block_write_gdal_image(outputName, left_trans_img, geo_opt);
#endif
      } //endif use_local_homography

      local_search_range = grow_bbox_to_int(local_search_range);
      // Expand local_search_range by 1. This is necessary since
      // m_sub_disp is integer-valued, and perhaps the search
      // range was supposed to be a fraction of integer bigger.
      local_search_range.expand(1);
      
      // Scale the search range to full-resolution
      local_search_range.min() = floor(elem_prod(local_search_range.min(),m_upscale_factor));
      local_search_range.max() = ceil (elem_prod(local_search_range.max(),m_upscale_factor));

      // If the user specified a search range limit, apply it here.
      if ((stereo_settings().search_range_limit.min() != Vector2i()) || 
          (stereo_settings().search_range_limit.max() != Vector2i())   ) {     
        local_search_range.crop(stereo_settings().search_range_limit);
        vw_out() << "\t--> Local search range constrained to: " << local_search_range << "\n";
      }

      VW_OUT(DebugMessage, "stereo") << "SeededCorrelatorView("
				     << bbox << ") local search range "
				     << local_search_range << " vs "
				     << stereo_settings().search_range << "\n";

    } else{ // seed mode == 0
      local_search_range = stereo_settings().search_range;
      VW_OUT(DebugMessage,"stereo") << "Searching with " << stereo_settings().search_range << "\n";
    }

    SemiGlobalMatcher::SgmSubpixelMode sgm_subpixel_mode = get_sgm_subpixel_mode();
    Vector2i sgm_search_buffer = stereo_settings().sgm_search_buffer;

    // Now we are ready to actually perform correlation
    const int rm_half_kernel = 5; // Filter kernel size used by CorrelationView
    if (use_local_homography){
		//<RM>: apply stereo to the aligned left and right tile
		typedef vw::stereo::PyramidCorrelationView<ImageView<InputPixelType>, ImageView<InputPixelType>, 
                                                   ImageView<vw::uint8     >, ImageView<vw::uint8     > > CorrView;
		CorrView corr_view( left_trans_img,   right_trans_img,
                          left_trans_mask,    right_trans_mask,
                          static_cast<vw::stereo::PrefilterModeType>(stereo_settings().pre_filter_mode),
                          stereo_settings().slogW,
                          fullres_hom != math::identity_matrix<3>() ? new_local_search_range : local_search_range, 
                          m_kernel_size,  m_cost_mode,
                          m_corr_timeout, m_seconds_per_op,
                          stereo_settings().xcorr_threshold,
                          stereo_settings().min_xcorr_level,
                          rm_half_kernel,
                          stereo_settings().corr_max_levels,
                          static_cast<vw::stereo::CorrelationAlgorithm>(stereo_settings().stereo_algorithm), 
                          stereo_settings().sgm_collar_size,
                          sgm_subpixel_mode, sgm_search_buffer, stereo_settings().corr_memory_limit_mb,
                          stereo_settings().corr_blob_filter_area,
                          stereo_settings().stereo_debug );
		ImageView<pixel_type> stereo_result = corr_view.prerasterize(bounding_box(left_trans_img));
#if DEBUG_RM
		cout << "[tile(" << H << "," << W << " Stereo done!" << endl;
#endif
      	ImageView<pixel_type> stereo_result_inv;
      	ImageView<vw::uint8> stereo_result_mask_inv;  
	  	ImageView<vw::uint8> stereo_result_mask = left_trans_mask;
		//<RM>: the disparity map needs to fit the original tile size therefore the transform needs to be undone
	  	ImageView< PixelMask<pixel_type> > stereo_result_masked_img_inv = transform (copy_mask(stereo_result.impl(), 																			  stereo_result_mask.impl()),
	               														  HomographyTransform(inverse(align_left_matrix)),
																		  newBBox.width(), newBBox.height());
      	stereo_result_inv  = apply_mask(stereo_result_masked_img_inv);
      	stereo_result_mask_inv = channel_cast_rescale<uint8>(select_channel(stereo_result_masked_img_inv, 2));
		//<RM>: remove the margin
		ImageView<pixel_type> stereo_result_corrected(bbox.width(), bbox.height());
		double marginMinX = bbox.min().x() == 0 ? 0 : margin;
		double marginMinY = bbox.min().y() == 0 ? 0 : margin;
		for(int j=0; j<bbox.height(); j++ ){
			for(int i=0; i<bbox.width(); i++ ){
				stereo_result_corrected(i,j)[0] = stereo_result_inv(i+marginMinX ,j+marginMinY)[0];
				stereo_result_corrected(i,j)[1] = stereo_result_inv(i+marginMinX, j+marginMinY)[1];
				if(stereo_result_mask_inv(i+marginMinX ,j+marginMinY))
					validate(stereo_result_corrected(i,j));
			}
		}
#if DEBUG_RM
		cout << "[tile(" << H << "," << W << " marginMinX = " << marginMinX << endl;
		cout << "[tile(" << H << "," << W << " marginMinY = " << marginMinY << endl;
		sprintf(outputName, "stereo_%d_%d.tif", H, W);
	 	block_write_gdal_image(outputName, stereo_result, geo_opt);
		sprintf(outputName, "stereoINV_%d_%d.tif", H, W);	
		block_write_gdal_image(outputName, stereo_result_inv, geo_opt);  
		sprintf(outputName, "stereoINVCorrected_%d_%d.tif", H, W);	
		block_write_gdal_image(outputName, stereo_result_corrected, geo_opt);  
		cout << "end of tile " << bbox << endl;	
#endif
	if(fullres_hom != math::identity_matrix<3>()) 
		cout << "local_search_range " << new_local_search_range << endl;
	else
		cout << "local_search_range " << local_search_range << endl;
		return prerasterize_type(stereo_result_corrected,-bbox.min().x(),-bbox.min().y(),cols(),rows() );
    }else{
      typedef vw::stereo::PyramidCorrelationView<ImageType, ImageType, MaskType, MaskType > CorrView;
      CorrView corr_view( m_left_image,   m_right_image,
                          m_left_mask,    m_right_mask,
                          static_cast<vw::stereo::PrefilterModeType>(stereo_settings().pre_filter_mode),
                          stereo_settings().slogW,
                          local_search_range,
                          m_kernel_size,  m_cost_mode,
                          m_corr_timeout, m_seconds_per_op,
                          stereo_settings().xcorr_threshold,
                          stereo_settings().min_xcorr_level,
                          rm_half_kernel,
                          stereo_settings().corr_max_levels,
                          static_cast<vw::stereo::CorrelationAlgorithm>(stereo_settings().stereo_algorithm), 
                          stereo_settings().sgm_collar_size,
                          sgm_subpixel_mode, sgm_search_buffer, stereo_settings().corr_memory_limit_mb,
                          stereo_settings().corr_blob_filter_area,
                          stereo_settings().stereo_debug );
	cout << "local_search_range " << local_search_range << endl;
      return corr_view.prerasterize(bbox);
    }
    
  } // End function prerasterize_helper

  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
}; // End class SeededCorrelatorView


/// Main stereo correlation function, called after parsing input arguments.
void stereo_correlation( ASPGlobalOptions& opt ) {

  // The first thing we will do is compute the low-resolution correlation.

  // Note that even when we are told to skip low-resolution correlation,
  // we must still go through the motions when seed_mode is 0, to be
  // able to get a search range, even though we don't write D_sub then.
  if (!stereo_settings().skip_low_res_disparity_comp || stereo_settings().seed_mode == 0)
    lowres_correlation(opt);

  if (stereo_settings().compute_low_res_disparity_only) 
    return; // Just computed the low-res disparity, so quit.

  vw_out() << "\n[ " << current_posix_time_string() << " ] : Stage 1 --> CORRELATION \n";

  read_search_range_from_dsub(opt);

  // If the user specified a search range limit, apply it here.
  if ((stereo_settings().search_range_limit.min() != Vector2i()) || 
      (stereo_settings().search_range_limit.max() != Vector2i())   ) {     
    stereo_settings().search_range.crop(stereo_settings().search_range_limit);
    vw_out() << "\t--> Detected search range constrained to: " << stereo_settings().search_range << "\n";
  }


  // Provide the user with some feedback of what we are actually going to use.
  vw_out()   << "\t--------------------------------------------------\n";
  vw_out()   << "\t   Kernel Size:    " << stereo_settings().corr_kernel << endl;
  if ( stereo_settings().seed_mode > 0 )
    vw_out() << "\t   Refined Search: " << stereo_settings().search_range << endl;
  else
    vw_out() << "\t   Search Range:   " << stereo_settings().search_range << endl;
  vw_out()   << "\t   Cost Mode:      " << stereo_settings().cost_mode << endl;
  vw_out(DebugMessage) << "\t   XCorr Threshold: " << stereo_settings().xcorr_threshold << endl;
  vw_out(DebugMessage) << "\t   Prefilter:       " << stereo_settings().pre_filter_mode << endl;
  vw_out(DebugMessage) << "\t   Prefilter Size:  " << stereo_settings().slogW << endl;
  vw_out() << "\t--------------------------------------------------\n";

  // Load up for the actual native resolution processing
  DiskImageView<PixelGray<float> > left_disk_image (opt.out_prefix+"-L.tif"),
                                   right_disk_image(opt.out_prefix+"-R.tif");
  DiskImageView<vw::uint8> Lmask(opt.out_prefix + "-lMask.tif"),
                           Rmask(opt.out_prefix + "-rMask.tif");
  ImageViewRef<PixelMask<Vector2f> > sub_disp;
  std::string dsub_file   = opt.out_prefix+"-D_sub.tif";
  std::string spread_file = opt.out_prefix+"-D_sub_spread.tif";
  
  if ( stereo_settings().seed_mode > 0 )
    sub_disp = DiskImageView<PixelMask<Vector2f> >(dsub_file);
  ImageViewRef<PixelMask<Vector2i> > sub_disp_spread;
  if ( stereo_settings().seed_mode == 2 ||  stereo_settings().seed_mode == 3 ){
    // D_sub_spread is mandatory for seed_mode 2 and 3.
    sub_disp_spread = DiskImageView<PixelMask<Vector2i> >(spread_file);
  }else if ( stereo_settings().seed_mode == 1 ){
    // D_sub_spread is optional for seed_mode 1, we use it only if it is provided.
    if (fs::exists(spread_file)) {
      try {
        sub_disp_spread = DiskImageView<PixelMask<Vector2i> >(spread_file);
      }
      catch (...) {}
    }
  }

  ImageView<Matrix3x3> local_hom;
  ImageView<float> verticalDisp;
  if ( stereo_settings().seed_mode > 0 && stereo_settings().use_local_homography ){
    string local_hom_file = opt.out_prefix + "-local_hom.txt";
    read_local_homographies(local_hom_file, local_hom);
  }
//<RM>: added support for local tranformations applied to both tiles
  ImageView<Matrix3x3> local_hom_L;
  ImageView<Matrix3x3> local_size; 
if ( stereo_settings().seed_mode > 0 && stereo_settings().use_local_homography ){
    string local_hom_file = opt.out_prefix + "-local_hom.txt";
//<RM>: TODO: just create empty files
    read_local_homographies(local_hom_file, local_hom_L);
	read_local_homographies(local_hom_file, local_size);
  }

  stereo::CostFunctionType cost_mode = get_cost_mode_value();
  Vector2i kernel_size    = stereo_settings().corr_kernel;
  BBox2i   trans_crop_win = stereo_settings().trans_crop_win;
  int      corr_timeout   = stereo_settings().corr_timeout;
  double   seconds_per_op = 0.0;
  if (corr_timeout > 0)
    seconds_per_op = calc_seconds_per_op(cost_mode, left_disk_image, right_disk_image, kernel_size);

  // TODO: Check!
  // Set up the reference to the stereo disparity code
  // - Processing is limited to trans_crop_win for use with parallel_stereo.
  ImageViewRef<PixelMask<Vector2f> > fullres_disparity =
    crop(SeededCorrelatorView( left_disk_image, right_disk_image, Lmask, Rmask,
                               sub_disp, sub_disp_spread, local_hom,
							   local_hom_L, local_size,  //Ricardo Monteiro
							   kernel_size, 
                               cost_mode, corr_timeout, seconds_per_op),  
         trans_crop_win);

  // With SGM, we must do the entire image chunk as one tile. Otherwise,
  // if it gets done in smaller tiles, there will be artifacts at tile boundaries.
  bool using_sgm = (stereo_settings().stereo_algorithm > vw::stereo::CORRELATION_WINDOW);
  /*if (using_sgm) {
    Vector2i image_size = bounding_box(fullres_disparity).size();
    int max_dim = std::max(image_size[0], image_size[1]);
    if (stereo_settings().corr_tile_size_ovr < max_dim)
      vw_throw(ArgumentErr()
               << "Error: SGM processing is not permitted with a tile size smaller than the image!\n"
               << "Value of --corr-tile-size is " << stereo_settings().corr_tile_size_ovr
               << " but image size is " << image_size << ".\n" 
               << "Increase --corr-tile-size so the entire image fits in one tile, or "
               << "use parallel_stereo. Not that making --corr-tile-size larger than 9000 or so may "
               << "cause GDAL to crash.\n\n");
  }*/
  
  switch(stereo_settings().pre_filter_mode){
  case 2:
    vw_out() << "\t--> Using LOG pre-processing filter with "
             << stereo_settings().slogW << " sigma blur.\n"; 
    break;
  case 1:
    vw_out() << "\t--> Using Subtracted Mean pre-processing filter with "
	           << stereo_settings().slogW << " sigma blur.\n";
    break;
  default:
    vw_out() << "\t--> Using NO pre-processing filter." << endl;
  }

  cartography::GeoReference left_georef;
  bool   has_left_georef = read_georeference(left_georef,  opt.out_prefix + "-L.tif");
  bool   has_nodata      = false;
  double nodata          = -32768.0;

  string d_file = opt.out_prefix + "-D.tif";
  vw_out() << "Writing: " << d_file << "\n";
  if (stereo_settings().stereo_algorithm > vw::stereo::CORRELATION_WINDOW) {
    // SGM performs subpixel correlation in this step, so write out floats.
    
    // Rasterize the image first as one block, then write it out using multiple blocks.
    // - If we don't do this, the output image file is not tiled and handles very slowly.
    // - This is possible because with SGM the image must be small enough to fit in memory.
//<RM>: TESTING
    //ImageView<PixelMask<Vector2f> > result = fullres_disparity;
    //opt.raster_tile_size = Vector2i(ASPGlobalOptions::rfne_tile_size(),ASPGlobalOptions::rfne_tile_size());
    //vw::cartography::block_write_gdal_image(d_file, result,
	//		        has_left_georef, left_georef,
	//		        has_nodata, nodata, opt,
	//		        TerminalProgressCallback("asp", "\t--> Correlation :") );
	vw::cartography::block_write_gdal_image(d_file, fullres_disparity,
       				has_left_georef, left_georef,
       				has_nodata, nodata, opt,
       				TerminalProgressCallback("asp", "\t--> Correlation :") );
			        
  } else {
    // Otherwise cast back to integer results to save on storage space.
    vw::cartography::block_write_gdal_image(d_file, 
              pixel_cast<PixelMask<Vector2i> >(fullres_disparity),
			        has_left_georef, left_georef,
			        has_nodata, nodata, opt,
			        TerminalProgressCallback("asp", "\t--> Correlation :") );
  }
//<RM>: overwrite transformations applied to right tile and also write to file the transformations applied to the left tile and the aligned tile size for each tile
if ( stereo_settings().seed_mode > 0 && stereo_settings().use_local_homography ){
    string local_hom_file = opt.out_prefix + "-local_hom.txt";
    write_local_homographies(local_hom_file, local_hom);
    string local_hom_L_file = opt.out_prefix + "-local_hom_L.txt";
	write_local_homographies(local_hom_L_file, local_hom_L);
	string local_size_file = opt.out_prefix + "-local_size.txt";
	write_local_homographies(local_size_file, local_size);
  }

  vw_out() << "\n[ " << current_posix_time_string() << " ] : CORRELATION FINISHED \n";

} // End function stereo_correlation

int main(int argc, char* argv[]) {

  //try { // TODO: Restore
    xercesc::XMLPlatformUtils::Initialize();

    stereo_register_sessions();

    bool verbose = false;
    vector<ASPGlobalOptions> opt_vec;
    string output_prefix;
    asp::parse_multiview(argc, argv, CorrelationDescription(),
			 verbose, output_prefix, opt_vec);
    ASPGlobalOptions opt = opt_vec[0];

    // TODO: Restore SGM restrictions!
    // TODO: Interate with stereo_blend/parallel_stereo

    // Leave the number of parallel block threads equal to the default unless we
    //  are using SGM in which case only one block at a time should be processed.
    // - Processing multiple blocks is possible, but it is better to use a larger blocks
    //   with more threads applied to the single block.
    // - Thread handling is still a little confusing because opt.num_threads is ONLY used
    //   to control the number of parallel image blocks written at a time.  Everything else
    //   reads directly from vw_settings().default_num_threads()
    const bool using_sgm = (stereo_settings().stereo_algorithm > vw::stereo::CORRELATION_WINDOW);
    opt.num_threads = vw_settings().default_num_threads();
    if (using_sgm)
      opt.num_threads = 1;

    // Integer correlator requires large tiles
    //---------------------------------------------------------
    int ts = stereo_settings().corr_tile_size_ovr;
    
    // GDAL block write sizes must be a multiple to 16 so if the input value is
    //  not a multiple of 16 increase it until it is.
    const int TILE_MULTIPLE = 16;
    if (ts % TILE_MULTIPLE != 0)
      ts = ((ts / TILE_MULTIPLE) + 1) * TILE_MULTIPLE;
      
    opt.raster_tile_size = Vector2i(ts, ts);

    // Internal Processes
    //---------------------------------------------------------
    stereo_correlation( opt );
  
    xercesc::XMLPlatformUtils::Terminate();
  //} ASP_STANDARD_CATCHES;

  return 0;
}

//<RM>: Added functions

//<RM>: piecewiseAlignment_affineepipolar - search for ip matches between the left and right tiles using homography_ip_matching1, apply piecewise affineepipolar alignment, check sanity of the tranformations and estimate new local search range
BBox2f piecewiseAlignment_affineepipolar( ImageView<float> tile_left_image,
										  ImageView<float> tile_right_image,	
										  BBox2i bbox,
										  Vector2i& left_size,
										  Vector2i& right_size,
										  vw::Matrix<double>& left_matrix,
										  vw::Matrix<double>& right_matrix,
										  BBox2f local_search_range)
{
	using namespace vw;
	//<RM>: piecewise alignment is only applied if it generates an average vertical disparity value of ..
	double threshPiecewiseAlignment = 2;
	//<RM>: average vertical disparity
	double avgDeltaY = -1.0;
	double threshRANSAC = 100.0; // 20 before
	//<RM>: multiplier for the estimated search range
	double threshSearchRange = 1;

	Matrix<double> H;
	double left_nodata_value  = numeric_limits<double>::quiet_NaN();
  	double right_nodata_value = numeric_limits<double>::quiet_NaN();
    std::vector<ip::InterestPoint> matched_ip1, matched_ip2;
	std::vector<ip::InterestPoint> matchedRANSAC_ip1,matchedRANSAC_ip2;
	std::vector<ip::InterestPoint> matchedRANSAC_final_ip1, matchedRANSAC_final_ip2;
	ip::InterestPoint aux_r_ip, aux_l_ip;
	char outputName[30];
	double newAvgDeltaY;
	double *newAvgDeltaYP = &newAvgDeltaY;
#if DEBUG_RM 
	int X = bbox.min().x()/ASPGlobalOptions::corr_tile_size();
	int Y = bbox.min().y()/ASPGlobalOptions::corr_tile_size();
	sprintf(outputName, "matches_%d_%d", Y, X);
#endif
    try {
		homography_ip_matching1( tile_left_image, tile_right_image,
                                           stereo_settings().ip_per_tile,
                                           outputName, threshRANSAC, 
                                           left_nodata_value, right_nodata_value,
					  					   matchedRANSAC_ip1, matchedRANSAC_ip2); 
	}catch(...){}
	//<RM>: estimate global alignment for this specific tile
	avgDeltaY = calcAverageDeltaY(matchedRANSAC_ip1, matchedRANSAC_ip2); 
#if DEBUG_RM 
	cout << "[tile(" << Y << "," << X << ") avgDeltaY after global alignment = " << avgDeltaY << "]" << endl;
#endif
	//<RM>: if the alignment can be improved
	if(avgDeltaY != -1 || avgDeltaY >= threshPiecewiseAlignment){
#if DEBUG_RM 
		ip::write_binary_match_file(outputName, matchedRANSAC_ip1, matchedRANSAC_ip2); 
		cout << "[tile(" << Y << "," << X << ")" << matchedRANSAC_ip1.size() << " matching points]" << endl;
#endif
		std::vector<Vector3> ransac_ip1 = iplist_to_vectorlist(matchedRANSAC_ip1), ransac_ip2 = iplist_to_vectorlist(matchedRANSAC_ip2);
		try {	
			left_size = affine_epipolar_rectification1(	left_size, right_size, matchedRANSAC_ip1, 															matchedRANSAC_ip2, left_matrix, right_matrix );
		} catch ( ... ) {
		  	left_matrix = math::identity_matrix<3>();
			right_matrix = math::identity_matrix<3>();
			return local_search_range;
		}
		//<RM>: check left_matrix and right_matrix
		if(!check_homography_matrix(left_matrix, right_matrix, ransac_ip1, ransac_ip2, avgDeltaY, bbox, newAvgDeltaYP)){
			left_matrix = math::identity_matrix<3>();
			right_matrix = math::identity_matrix<3>();
			return local_search_range;
		}
	//<RM>: if the alignment cannot be improved
	}else{ 
		left_matrix = math::identity_matrix<3>();
		right_matrix = math::identity_matrix<3>();
		return local_search_range;
	}
	//<RM>: in case the affine_epipolar_rectification1 returned identity matrix
	if(left_matrix == math::identity_matrix<3>() && right_matrix == math::identity_matrix<3>())
		return local_search_range;

	//<RM>: estimate new search range
	//return calcSearchRange(matchedRANSAC_ip1, matchedRANSAC_ip2, left_matrix, right_matrix, threshSearchRange);
	return calcSearchRange(matchedRANSAC_ip1, matchedRANSAC_ip2, left_matrix, right_matrix, newAvgDeltaY);
}

//<RM>: piecewiseAlignment_homography - Alternative method to Affine Epipolar - Not very good
BBox2f piecewiseAlignment_homography(ImageView<float> tile_left_image,
									 ImageView<float> tile_right_image,	
									 BBox2i bbox,
									 Vector2i& left_size,
									 Vector2i& right_size,
									 vw::Matrix<double>& left_matrix,
									 vw::Matrix<double>& right_matrix,
									 BBox2f local_search_range)
{
	using namespace vw;
	//<RM>: piecewise alignment is only applied if it generates an average vertical disparity value of ..
	double threshPiecewiseAlignment = 2;
	//<RM>: average vertical disparity
	double avgDeltaY = -1.0;
	double threshRANSAC = 100.0; // 20 before
	//<RM>: multiplier for the estimated search range
	double threshSearchRange = 1;

	Matrix<double> H;
	double left_nodata_value  = numeric_limits<double>::quiet_NaN();
  	double right_nodata_value = numeric_limits<double>::quiet_NaN();
    std::vector<ip::InterestPoint> matched_ip1, matched_ip2;
	std::vector<ip::InterestPoint> matchedRANSAC_ip1,matchedRANSAC_ip2;
	std::vector<ip::InterestPoint> matchedRANSAC_final_ip1, matchedRANSAC_final_ip2;
	ip::InterestPoint aux_r_ip, aux_l_ip;
	char outputName[30];
	double newAvgDeltaY;
	double *newAvgDeltaYP = &newAvgDeltaY;
#if DEBUG_RM 
	int X = bbox.min().x()/ASPGlobalOptions::corr_tile_size();
	int Y = bbox.min().y()/ASPGlobalOptions::corr_tile_size();
	sprintf(outputName, "matches_%d_%d", Y, X);
#endif
    try {
		homography_ip_matching1( tile_left_image, tile_right_image,
                                           stereo_settings().ip_per_tile,
                                           outputName, threshRANSAC, 
                                           left_nodata_value, right_nodata_value,
					  					   matchedRANSAC_ip1, matchedRANSAC_ip2); 
	}catch(...){}
	//<RM>: estimate global alignment for this specific tile
	avgDeltaY = calcAverageDeltaY(matchedRANSAC_ip1, matchedRANSAC_ip2); 
#if DEBUG_RM 
	cout << "[tile(" << Y << "," << X << ") avgDeltaY after global alignment = " << avgDeltaY << "]" << endl;
#endif
	//<RM>: if the alignment can be improved
	if(avgDeltaY != -1 || avgDeltaY >= threshPiecewiseAlignment){
#if DEBUG_RM 
		ip::write_binary_match_file(outputName, matchedRANSAC_ip1, matchedRANSAC_ip2); 
		cout << "[tile(" << Y << "," << X << ")" << matchedRANSAC_ip1.size() << " matching points]" << endl;
#endif
		std::vector<Vector3> ransac_ip1 = iplist_to_vectorlist(matchedRANSAC_ip1), ransac_ip2 = iplist_to_vectorlist(matchedRANSAC_ip2);
		try {	
 			left_size = homography_rectification1( false, tile_left_image.get_size(), tile_right_image.get_size(),
                                	   				matchedRANSAC_ip1, matchedRANSAC_ip2, left_matrix, right_matrix, 
									   				threshRANSAC, avgDeltaY, bbox );
		} catch ( ... ) {
		  	left_matrix = math::identity_matrix<3>();
			right_matrix = math::identity_matrix<3>();
			return local_search_range;
		}
		//<RM>: check left_matrix and right_matrix
		if(!check_homography_matrix(left_matrix, right_matrix, ransac_ip1, ransac_ip2, avgDeltaY, bbox, newAvgDeltaYP)){
			left_matrix = math::identity_matrix<3>();
			right_matrix = math::identity_matrix<3>();
			return local_search_range;
		}
	//<RM>: if the alignment cannot be improved
	}else{ 
		left_matrix = math::identity_matrix<3>();
		right_matrix = math::identity_matrix<3>();
		return local_search_range;
	}
	//<RM>: in case the affine_epipolar_rectification1 returned identity matrix
	if(left_matrix == math::identity_matrix<3>() && right_matrix == math::identity_matrix<3>())
		return local_search_range;

	//<RM>: estimate new search range
	//return calcSearchRange(matchedRANSAC_ip1, matchedRANSAC_ip2, left_matrix, right_matrix, threshSearchRange);
	return calcSearchRange(matchedRANSAC_ip1, matchedRANSAC_ip2, left_matrix, right_matrix, newAvgDeltaY);
}

//<RM>: affine_epipolar_rectification1 (copy of the original function) 
Vector2i affine_epipolar_rectification1( Vector2i const& left_size,
                                 Vector2i const& right_size,
                                 std::vector<ip::InterestPoint> const& ip1,
                                 std::vector<ip::InterestPoint> const& ip2,
                                 Matrix<double>& left_matrix,
                                 Matrix<double>& right_matrix ) {
    // Create affine fundamental matrix
    Matrix<double> fund = linear_affine_fundamental_matrix( ip1, ip2 );

    // Solve for rotation matrices
    double Hl = sqrt( fund(2,0)*fund(2,0) + fund(2,1)*fund(2,1) );
    double Hr = sqrt( fund(0,2)*fund(0,2) + fund(1,2)*fund(1,2) );
    Vector2 epipole(-fund(2,1),fund(2,0)), epipole_prime(-fund(1,2),fund(0,2));
    if ( epipole.x() < 0 )
      epipole = -epipole;
    if ( epipole_prime.x() < 0 )
      epipole_prime = -epipole_prime;
    epipole.y() = -epipole.y();
    epipole_prime.y() = -epipole_prime.y();

    left_matrix = math::identity_matrix<3>();
    right_matrix = math::identity_matrix<3>();
    left_matrix(0,0) = epipole[0]/Hl;
    left_matrix(0,1) = -epipole[1]/Hl;
    left_matrix(1,0) = epipole[1]/Hl;
    left_matrix(1,1) = epipole[0]/Hl;
    right_matrix(0,0) = epipole_prime[0]/Hr;
    right_matrix(0,1) = -epipole_prime[1]/Hr;
    right_matrix(1,0) = epipole_prime[1]/Hr;
    right_matrix(1,1) = epipole_prime[0]/Hr;

    // Solve for ideal scaling and translation
    solve_y_scaling( ip1, ip2, left_matrix, right_matrix );

    // Solve for ideal shear, scale, and translation of X axis
    solve_x_shear( ip1, ip2, left_matrix, right_matrix );

    // Work out the ideal render size.
    BBox2i output_bbox, right_bbox;
    output_bbox.grow( subvector(left_matrix*Vector3(0,0,1),0,2) );
    output_bbox.grow( subvector(left_matrix*Vector3(left_size.x(),0,1),0,2) );
    output_bbox.grow( subvector(left_matrix*Vector3(left_size.x(),left_size.y(),1),0,2) );
    output_bbox.grow( subvector(left_matrix*Vector3(0,left_size.y(),1),0,2) );
    right_bbox.grow( subvector(right_matrix*Vector3(0,0,1),0,2) );
    right_bbox.grow( subvector(right_matrix*Vector3(right_size.x(),0,1),0,2) );
    right_bbox.grow( subvector(right_matrix*Vector3(right_size.x(),right_size.y(),1),0,2) );
    right_bbox.grow( subvector(right_matrix*Vector3(0,right_size.y(),1),0,2) );
  //  output_bbox.crop( right_bbox );
	output_bbox.grow( right_bbox ); //<RM>: To avoid part of the tile being cropped

    left_matrix(0,2) -= output_bbox.min().x();
    right_matrix(0,2) -= output_bbox.min().x();
    left_matrix(1,2) -= output_bbox.min().y();
    right_matrix(1,2) -= output_bbox.min().y();

	//<RM>: Absolute limit for the output tile size (testing with a tile size of 1024 * 1.8)
	if(output_bbox.width() > 10000 || output_bbox.height() > 10000)
	{
		left_matrix = math::identity_matrix<3>();
		right_matrix = math::identity_matrix<3>();
		return left_size;
	}

    return Vector2i( output_bbox.width(), output_bbox.height() );
  }

//<RM>: homography_ip_matching1 (copy of the original function) - the function was modified to be used for piecewise alignment 
  template <class Image1T, class Image2T>
  bool homography_ip_matching1( vw::ImageViewBase<Image1T> const& image1,
			       vw::ImageViewBase<Image2T> const& image2,
			       int ip_per_tile,
			       std::string const& output_name,
			       int inlier_threshold,
			       double nodata1,
			       double nodata2,
			       std::vector<ip::InterestPoint>& final_ip1,
 			       std::vector<ip::InterestPoint>& final_ip2) {

    using namespace vw;

    std::vector<ip::InterestPoint> matched_ip1, matched_ip2;
    detect_match_ip( matched_ip1, matched_ip2,
		     image1.impl(), image2.impl(),
		     ip_per_tile,
		     nodata1, nodata2 );

    if ( matched_ip1.size() == 0 || matched_ip2.size() == 0 )
      return false;
    std::vector<Vector3> ransac_ip1 = iplist_to_vectorlist(matched_ip1),
			 ransac_ip2 = iplist_to_vectorlist(matched_ip2);
    std::vector<size_t> indices;
    try {
      typedef math::RandomSampleConsensus<math::HomographyFittingFunctor, math::InterestPointErrorMetric> RansacT;
      const int    MIN_NUM_OUTPUT_INLIERS = ransac_ip1.size()/2;
      const int    NUM_ITERATIONS         = 100;
      RansacT ransac( math::HomographyFittingFunctor(),
		      math::InterestPointErrorMetric(), NUM_ITERATIONS,
		      inlier_threshold,
		      MIN_NUM_OUTPUT_INLIERS, true
		      );
      Matrix<double> H(ransac(ransac_ip2,ransac_ip1)); // 2 then 1 is used here for legacy reasons
      indices = ransac.inlier_indices(H,ransac_ip2,ransac_ip1);
    } catch (const math::RANSACErr& e ) {
      return false;
    }

    BOOST_FOREACH( size_t& index, indices ) {
      final_ip1.push_back(matched_ip1[index]);
      final_ip2.push_back(matched_ip2[index]);
    }
#if DEBUG_RM
    ip::write_binary_match_file(output_name, final_ip1, final_ip2);
#endif
    return true;
  }

//<RM>: homography_rectification1 (copy of the original function) - the function was modified to be used for piecewise alignment 
Vector2i homography_rectification1( bool adjust_left_image_size,
			    					Vector2i const& left_size,
			    					Vector2i const& right_size,
			    					std::vector<ip::InterestPoint> const& left_ip,
			    					std::vector<ip::InterestPoint> const& right_ip,
			    					vw::Matrix<double>& left_matrix,
			    					vw::Matrix<double>& right_matrix,
			    					double threshRANSAC,
			    					double minAvgDeltaY, 
									BBox2i bbox ) {
    // Reformat the interest points for RANSAC
    std::vector<Vector3>  right_copy = iplist_to_vectorlist(right_ip),
			  left_copy  = iplist_to_vectorlist(left_ip);
    
    // Use RANSAC to determine a good homography transform between the images
    math::RandomSampleConsensus<math::HomographyFittingFunctor, math::InterestPointErrorMetric>
      ransac( math::HomographyFittingFunctor(),
	      math::InterestPointErrorMetric(),
	      100, // num iter
	      threshRANSAC,
	      left_copy.size()/2 // min output inliers
	      );
	
    Matrix<double> H = ransac(right_copy, left_copy);
    std::vector<size_t> indices = ransac.inlier_indices(H, right_copy, left_copy);

    if(check_homography_matrix(H, left_copy, right_copy, indices, minAvgDeltaY, bbox)){
    // Set right to a homography that has been refined just to our inliers
    	left_matrix  = math::identity_matrix<3>();
    	right_matrix = math::HomographyFittingFunctor()(right_copy, left_copy, H);
    }else{
		left_matrix  = math::identity_matrix<3>();
    	right_matrix = math::identity_matrix<3>();
    }

// Ricardo Monteiro // calculate average deltat Y after piecewise alignment
    std::vector<ip::InterestPoint> left_points, right_points;
    int j = 0;
    for( size_t i = 0; i < left_ip.size(); i++ ){
	if(i == indices[j]){
	    left_points.push_back(left_ip[i]);
	    right_points.push_back(right_ip[i]);
	    j++;
	}
    }
    double averageDeltaY = calcAverageDeltaY(left_points, right_points);
    cout << "[" << averageDeltaY << " averageDeltaY after piecewise alignment]" << endl;

    // Work out the ideal render size
    BBox2i output_bbox, right_bbox;
    output_bbox.grow( Vector2i(0,0) );
    output_bbox.grow( Vector2i(left_size.x(),0) );
    output_bbox.grow( Vector2i(0,left_size.y()) );
    output_bbox.grow( left_size );

    if (adjust_left_image_size){
      // Crop the left and right images to the shared region. This is
      // done for efficiency.  It may not be always desirable though,
      // as in this case we lose the one-to-one correspondence between
      // original input left image pixels and output disparity/point
      // cloud pixels.
      Vector3 temp = right_matrix*Vector3(0,0,1);
      temp /= temp.z();
      right_bbox.grow( subvector(temp,0,2) );
      temp = right_matrix*Vector3(right_size.x(),0,1);
      temp /= temp.z();
      right_bbox.grow( subvector(temp,0,2) );
      temp = right_matrix*Vector3(0,right_size.y(),1);
      temp /= temp.z();
      right_bbox.grow( subvector(temp,0,2) );
      temp = right_matrix*Vector3(right_size.x(),right_size.y(),1);
      temp /= temp.z();
      right_bbox.grow( subvector(temp,0,2) );

      output_bbox.crop( right_bbox );

      //  Move the ideal render size to be aligned up with origin
      left_matrix (0,2) -= output_bbox.min().x();
      right_matrix(0,2) -= output_bbox.min().x();
      left_matrix (1,2) -= output_bbox.min().y();
      right_matrix(1,2) -= output_bbox.min().y();
    }

    return Vector2i( output_bbox.width(), output_bbox.height() );
  }

//<RM>: check_homography_matrix - sanity check for both left and right matrices (NOT being used)
bool check_homography_matrix(Matrix<double>       const& H,
                               std::vector<Vector3> const& left_points,
                               std::vector<Vector3> const& right_points,
                               std::vector<size_t>  const& indices,
                               double minAvgDeltaY, BBox2i bbox
                               ){

    // Sanity checks. If these fail, most likely the two images are too different
    // for stereo to succeed.
    /*if ( indices.size() < std::min( right_points.size(), left_points.size() )/2 ){
      vw_out(WarningMessage) << "InterestPointMatching: The number of inliers is less "
                             << "than 1/2 of the number of points. The inputs may be invalid.\n";
        return false;
    }*/

    double det = fabs(H(0, 0)*H(1, 1) - H(0, 1)*H(1, 0));
    if (det <= 0.5 || det >= 2.0){
      vw_out(WarningMessage) << "InterestPointMatching: The determinant of the 2x2 submatrix "
                             << "of the homography matrix " << H << " is " << det
                             << ". There could be a large scale discrepancy among the input images "
                             << "or the inputs may be an invalid stereo pair.\n";
        return false;
    }

    // check if the avgDeltaY after piecewise alignment is better than the minAvgDeltaY
    std::vector<Vector3> right_ip;
    std::vector<Vector3> left_ip;
        std::vector<ip::InterestPoint> r_ip, l_ip;
        ip::InterestPoint aux_r_ip, aux_l_ip;
    double avgDeltaY = -1;
    int ts = ASPGlobalOptions::corr_tile_size();
    for(size_t i = 0; i < right_points.size(); i++)
    {
                //cout << " ip matchings " << right_points[i].y() << " " << right_points[i].x() << "]" << endl; // DEBUG_RM
        right_ip.push_back(H * Vector3(right_points[i].x(), right_points[i].y(), 1));
        left_ip.push_back(Vector3(left_points[i].x(), left_points[i].y(), 1));
        // Normalize the coordinates, but don't divide by 0
        if (right_ip[i].z() == 0 || left_ip[i].z() == 0)
            continue;
        right_ip[i] /= right_ip[i].z();
        left_ip[i] /= left_ip[i].z();
#if DEBUG_RM
        aux_l_ip.x = left_ip[i].x() - bbox.min().x(); 
        aux_l_ip.y = left_ip[i].y() - bbox.min().y(); 
        aux_r_ip.x = right_ip[i].x() - bbox.min().x(); 
        aux_r_ip.y = right_ip[i].y() - bbox.min().y(); 
        r_ip.push_back(aux_r_ip); 
        l_ip.push_back(aux_l_ip); 
        //cout << " ip matchings after H " << right_ip[i].y() << " " << right_ip[i].x() << "]" << endl; // DEBUG_RM
#endif
    }
    avgDeltaY = calcAverageDeltaY(left_ip, right_ip);
    cout << "[tile(" << bbox.min().y()/ts << "," << bbox.min().x()/ts << ") avgDeltaY after piecewise alignment = " << avgDeltaY << "]" << endl;
#if DEBUG_RM
        char outputName[30];
        int X = bbox.min().x()/ASPGlobalOptions::corr_tile_size(); 
        int Y = bbox.min().y()/ASPGlobalOptions::corr_tile_size(); 
        sprintf(outputName, "matches_after_H_%d_%d", Y, X);
        ip::write_binary_match_file(outputName, l_ip, r_ip); 
#endif
    if(avgDeltaY == -1 || avgDeltaY >= minAvgDeltaY)
        return false;

    return true;

  }

//<RM>: check_homography_matrix - sanity check for both left and right matrices 
bool check_homography_matrix(	Matrix<double>       const& left_matrix,
								Matrix<double>       const& right_matrix,
			       				std::vector<Vector3> const& left_points,
			       				std::vector<Vector3> const& right_points,
			       				double minAvgDeltaY, 
								BBox2i bbox, double *avgDeltaY
								){
    double det = fabs(left_matrix(0, 0)*left_matrix(1, 1) - left_matrix(0, 1)*left_matrix(1, 0));
    if (det <= 0.5 || det >= 2.0){
      vw_out(WarningMessage) << "InterestPointMatching: The determinant of the 2x2 submatrix "
                             << "of the homography matrix " << left_matrix << " is " << det
                             << ". There could be a large scale discrepancy among the input images "
                             << "or the inputs may be an invalid stereo pair.\n";
	return false;
    }
	det = fabs(right_matrix(0, 0)*right_matrix(1, 1) - right_matrix(0, 1)*right_matrix(1, 0));
    if (det <= 0.5 || det >= 2.0){
      vw_out(WarningMessage) << "InterestPointMatching: The determinant of the 2x2 submatrix "
                             << "of the homography matrix " << right_matrix << " is " << det
                             << ". There could be a large scale discrepancy among the input images "
                             << "or the inputs may be an invalid stereo pair.\n";
	return false;
    }

    //<RM>: check if the avgDeltaY after piecewise alignment is better than the minAvgDeltaY
    std::vector<Vector3> right_ip; 
    std::vector<Vector3> left_ip;
	std::vector<ip::InterestPoint> r_ip, l_ip;
	ip::InterestPoint aux_r_ip, aux_l_ip;
    //double avgDeltaY = -1;
    for(size_t i = 0; i < right_points.size(); i++)
    { 
        right_ip.push_back(right_matrix * Vector3(right_points[i].x(), right_points[i].y(), 1));
        left_ip.push_back(left_matrix * Vector3(left_points[i].x(), left_points[i].y(), 1));
		//<RM>: Normalize the coordinates, but don't divide by 0
        if (right_ip[i].z() == 0 || left_ip[i].z() == 0) 
            continue;
        right_ip[i] /= right_ip[i].z();
        left_ip[i] /= left_ip[i].z();
#if DEBUG_RM
		aux_l_ip.x = left_ip[i].x(); 
    	aux_l_ip.y = left_ip[i].y() ; 
    	aux_r_ip.x = right_ip[i].x(); 
    	aux_r_ip.y = right_ip[i].y(); 
		r_ip.push_back(aux_r_ip); 
		l_ip.push_back(aux_l_ip);
#endif
    }
	//<RM>: calculate average vertical disparity after piecewise alignment
    *avgDeltaY = calcAverageDeltaY(left_ip, right_ip);
#if DEBUG_RM
	int ts = ASPGlobalOptions::corr_tile_size();
	char outputName[30]; // DEBUG_RM
	int X = bbox.min().x()/ASPGlobalOptions::corr_tile_size(); 
	int Y = bbox.min().y()/ASPGlobalOptions::corr_tile_size(); 
	sprintf(outputName, "matches_after_H_%d_%d", Y, X); 
	ip::write_binary_match_file(outputName, l_ip, r_ip);
	cout << "[tile(" << bbox.min().y()/ts << "," << bbox.min().x()/ts << ") avgDeltaY after piecewise alignment = " << *avgDeltaY << "]" << endl;
#endif
    if(*avgDeltaY == -1 || *avgDeltaY >= minAvgDeltaY)
        return false;

    return true;

  }

//<RM>: calcAverageDeltaY (vector<ip::InterestPoint>) - estimate average alignment based on ip matches
double calcAverageDeltaY(std::vector<ip::InterestPoint> const& left_points, std::vector<ip::InterestPoint> const& right_points)
{
    double accuDiff = 0;

    if(left_points.size()){
        for ( size_t i = 0; i < left_points.size(); i++ )
	    accuDiff += abs(left_points[i].y - right_points[i].y);
        return accuDiff/left_points.size(); //<RM>: average
    }else
	return -1; //<RM>: not valid
}
//<RM>: calcAverageDeltaY (Vector3) - estimate average alignment based on ip matches
double calcAverageDeltaY(std::vector<Vector3> const& left_points, std::vector<Vector3> const& right_points)
{
    double accuDiff = 0;

    if(left_points.size()){
        for ( size_t i = 0; i < left_points.size(); i++ )
	    accuDiff += abs(left_points[i].y() - right_points[i].y());
        return accuDiff/left_points.size(); //<RM>: average
    }else
	return -1; //<RM>: not valid
}

//<RM>: calcSearchRange - estimate search range and multiply it by threshold 
BBox2f calcSearchRange(std::vector<ip::InterestPoint> const& left_ip, std::vector<ip::InterestPoint> const& right_ip, Matrix<double> const& left_matrix, Matrix<double> const& right_matrix, double multi)
{
	std::vector<int> diffY, diffX;
	int maxDiffY = 0, maxDiffX = 0, minDiffY = 0, minDiffX = 0; 
	std::vector<Vector3> trans_left_points, trans_right_points;
	int horizontalTH = 2;

	//<RM>: transform ip matches
	for(size_t i = 0; i < left_ip.size(); i++)  
    { 
        trans_right_points.push_back(right_matrix * Vector3(right_ip[i].x, right_ip[i].y, 1));
        trans_left_points.push_back(left_matrix * Vector3(left_ip[i].x, left_ip[i].y, 1));
        if (trans_right_points[i].z() == 0 || trans_left_points[i].z() == 0) 
            continue;
        trans_right_points[i] /= trans_right_points[i].z();
        trans_left_points[i] /= trans_left_points[i].z();
    }
	//<RM>: gen list of diff (right - left)
	for ( size_t i = 0; i < trans_right_points.size(); i++ ){ 
		diffY.push_back(trans_right_points[i].y() - trans_left_points[i].y());
		diffX.push_back(trans_right_points[i].x() - trans_left_points[i].x());
	}
	//<RM>: get min and max
	for ( size_t i = 0; i < trans_right_points.size(); i++ ){
		if(diffY[i] < minDiffY)
			minDiffY = diffY[i];
		if(diffX[i] < minDiffX)
			minDiffX = diffX[i];
		if(diffY[i] > maxDiffY)
			maxDiffY = diffY[i];
		if(diffX[i] > maxDiffX)
			maxDiffX = diffX[i];
	}
	cout << "avgDeltaY: " << multi <<  endl;
	//cout << "floor(-2*multi): " << floor(-2*multi) <<  endl;
	//cout << "ceil(4*multi+1): " << ceil(4*multi+1) <<  endl;
	return BBox2f(horizontalTH*minDiffX, floor(-2*multi), (horizontalTH*maxDiffX) - (2*minDiffX), ceil(4*multi+1));
	//return BBox2f(multi * minDiffX, multi * minDiffY, (multi * maxDiffX) - (multi * minDiffX), (multi * maxDiffY) - (multi * minDiffY));
	//return BBox2f(2 * minDiffX, 2 * minDiffY, (2 * maxDiffX) - (2 * minDiffX), (2 * maxDiffY) - (2 * minDiffY));
	//return BBox2f(multi * minDiffX, 0, (multi * maxDiffX) - (multi * minDiffX), 1);

}

