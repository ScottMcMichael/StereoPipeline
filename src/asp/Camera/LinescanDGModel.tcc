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


template <class PositionFuncT, class VelocityFuncT, class PoseFuncT, class TimeFuncT>
typename LinescanDGModel<PositionFuncT, VelocityFuncT, PoseFuncT, TimeFuncT>::LinescanLMA::result_type
LinescanDGModel<PositionFuncT, VelocityFuncT, PoseFuncT, TimeFuncT>::LinescanLMA
::operator()( domain_type const& y ) const {
  double t = m_model->m_time_func( y[0] );

  // Rotate the point into our camera's frame
  vw::Vector3 pt = inverse( m_model->m_pose_func(t) ).rotate( m_point - m_model->m_position_func(t) );
  pt *= m_model->m_focal_length / pt.z(); // Rescale to pixel units
  result_type result(1);
  result[0] = pt.y() -
    m_model->m_detector_origin[1]; // Error against the location
				   // of the detector
  return result;
}

template <class PositionFuncT, class VelocityFuncT, class PoseFuncT, class TimeFuncT>
vw::Vector2 LinescanDGModel<PositionFuncT, VelocityFuncT, PoseFuncT, TimeFuncT>
::point_to_pixel(vw::Vector3 const& point) const {
  return point_to_pixel(point, -1);
}

// Here we use an initial guess for the line number
template <class PositionFuncT, class VelocityFuncT, class PoseFuncT, class TimeFuncT>
vw::Vector2 LinescanDGModel<PositionFuncT, VelocityFuncT, PoseFuncT, TimeFuncT>
::point_to_pixel(vw::Vector3 const& point, double starty) const {

  if (!m_correct_velocity_aberration)
    return point_to_pixel_uncorrected(point, starty);
  return point_to_pixel_corrected(point, starty);
}

namespace linescan {
  const double ABS_TOL = 1e-16;
  const double REL_TOL = 1e-16;
  const int    MAX_ITERATIONS = 1e+5;
}

template <class PositionFuncT, class VelocityFuncT, class PoseFuncT, class TimeFuncT>
vw::Vector2 LinescanDGModel<PositionFuncT, VelocityFuncT, PoseFuncT, TimeFuncT>
::point_to_pixel_uncorrected(vw::Vector3 const& point, double starty) const {

  using namespace vw;

  // Solve for the correct line number to use
  LinescanLMA model( this, point );
  int status;
  Vector<double> objective(1), start(1);
  start[0] = m_image_size.y()/2;

  // Use a refined guess, if available
  if (starty >= 0)
    start[0] = starty;

  Vector<double> solution
    = math::levenberg_marquardt(model, start, objective, status,
                                linescan::ABS_TOL,
                                linescan::REL_TOL,
                                linescan::MAX_ITERATIONS);

  VW_ASSERT( status > 0,
	     camera::PointToPixelErr() << "Unable to project point into LinescanDG model" );

  // Solve for sample location
  double  t  = m_time_func( solution[0] );
  Vector3 pt = inverse( m_pose_func(t) ).rotate( point - m_position_func(t) );
  pt *= m_focal_length / pt.z();

  return vw::Vector2(pt.x() - m_detector_origin[0], solution[0]);
}

template <class PositionFuncT, class VelocityFuncT, class PoseFuncT, class TimeFuncT>
vw::Vector2 LinescanDGModel<PositionFuncT, VelocityFuncT, PoseFuncT, TimeFuncT>
::point_to_pixel_corrected(vw::Vector3 const& point, double starty) const {

  using namespace vw;

  // Use the uncorrected function call to get an initial estimate
  Vector2 start = point_to_pixel_uncorrected(point, starty);

  // Improve the estimate using a model including the correction
  LinescanCorrLMA model( this, point );
  int status;
  Vector3 objective(0, 0, 0);
  Vector2 solution
    = math::levenberg_marquardt(model, start, objective, status,
                                linescan::ABS_TOL,
                                linescan::REL_TOL,
                                linescan::MAX_ITERATIONS);
  VW_ASSERT( status > 0,
	     camera::PointToPixelErr() << "Unable to project point into LinescanDG model" );

  return solution;
}

template <class PositionFuncT, class VelocityFuncT, class PoseFuncT, class TimeFuncT>
vw::Vector3 LinescanDGModel<PositionFuncT, VelocityFuncT, PoseFuncT, TimeFuncT>
::pixel_to_vector(vw::Vector2 const& pix) const {

  using namespace vw;

  // Compute local vector from the pixel out of the sensor
  // - m_detector_origin and m_focal_length have been converted into units of pixels
  Vector3 local_vec(pix[0]+m_detector_origin[0], m_detector_origin[1], m_focal_length);
  // Put the local vector in world coordinates using the pose information.
  Vector3 pix_to_vec = normalize(camera_pose(pix).rotate(local_vec));

  if (!m_correct_velocity_aberration) return pix_to_vec;


  // Correct for atmospheric refraction
  // - From "Atmospheric Correction and Velocity Aberration for Physical Sensor 
  //         Modeling of High-Resolution Satellite Images" by Oh and Lee

  // Get some information
  Vector3 cam_ctr              = camera_center(pix); // ECF camera coords
  Vector3 cam_vel              = camera_velocity(pix); // The camera's velocity in its own frame
  double  earth_ctr_to_cam     = norm_2(cam_ctr);    // Distance in meters from cam to earth center
  double  earth_rad            = 6371000.0; // TODO: Vary by location?
  double  cam_to_earth_surface = earth_ctr_to_cam - earth_rad;

    
  // Compute atmospheric condition constant K
  // TODO: Double check orig paper
  double H      = cam_to_earth_surface / 1000.0; // In kilometers
  double h      = 0; // Height of target point over surface // TODO: Use this! (in kilometers)
  double h_diff = H - h;
  double p1     = (2335.0 / h_diff)*pow(1 - 0.02257*h, 5.256);
  double p2     = pow(0.8540, H-11.0) * (82.2 - 521.0/h_diff);
  double K      = (p1 - p2) * 10e-6;
  
  // Compute angle alpha and correction angle
  // TODO: Should be a function call
  Vector3 cam_ctr_norm = normalize(cam_ctr);
  Vector3 u            = pix_to_vec;
  double  alpha        = acos(dot_prod(cam_ctr_norm,u) / (norm_2(cam_ctr_norm)*norm_2(u)) );  
  double  delta_alpha  = K * tan(alpha);

  
  // Rotate the vector by delta_alpha
  Vector3 rotation_axis = normalize(cross_prod(u, cam_ctr_norm));
  Quaternion<double> refraction_rotation(rotation_axis, delta_alpha);
  Vector3 u_prime = refraction_rotation.rotate(u);
  double  check_delta  = acos(dot_prod(u_prime,u) / (norm_2(u_prime)*norm_2(u)) );
/*
  vw_out() << "Computed K = " << K << std::endl;
  vw_out() << "Computed rot axis = " << rotation_axis << std::endl;
  vw_out() << "Computed delta_alpha = " << delta_alpha << std::endl;
  vw_out() << "Computed check = " << check_delta << std::endl;  
  vw_out() << "Computed u = " << u << std::endl;
  vw_out() << "Computed u prime = " << u_prime << std::endl;
*/
  pix_to_vec = u_prime;

  // Correct for velocity aberration
  // - If this is not done the positions will be well off the corresponding RPC locations.
  // - TODO: This should be an option for other satellite cameras, it is not
  //         specific to DG.

  // 1. Find the distance from the camera to the first
  // intersection of the current ray with the Earth surface.
  double  cam_angle_cos    = dot_prod(pix_to_vec, -normalize(cam_ctr));
  double  len_cos          = earth_ctr_to_cam*cam_angle_cos;
  double  cam_to_surface   = len_cos - sqrt(earth_rad*earth_rad
                              					    + len_cos*len_cos
                              					    - earth_ctr_to_cam*earth_ctr_to_cam);

  // 2. Correct the camera velocity due to the fact that the Earth
  // rotates around its axis.
  double seconds_in_day = 86164.0905;
  Vector3 earth_rotation_vec(0.0, 0.0, 2*M_PI/seconds_in_day); // Radians per second
  Vector3 cam_vel_corr1 = cam_vel - cam_to_surface * cross_prod(earth_rotation_vec, pix_to_vec);

  // 3. Find the component of the camera velocity orthogonal to the
  // direction the camera is pointing to.
  Vector3 cam_vel_corr2 = cam_vel_corr1 - dot_prod(cam_vel_corr1, pix_to_vec) * pix_to_vec;

  // 4. Correct direction for velocity aberration due to the speed of light.
  double light_speed = 299792458.0;
  Vector3 corr_pix_to_vec = pix_to_vec - cam_vel_corr2/light_speed;
  return normalize(corr_pix_to_vec);
}

template <class PositionFuncT, class VelocityFuncT, class PoseFuncT, class TimeFuncT>
vw::camera::PinholeModel LinescanDGModel<PositionFuncT, VelocityFuncT, PoseFuncT, TimeFuncT>
::linescan_to_pinhole(double y) const{

  double t = m_time_func( y );
  return vw::camera::PinholeModel(m_position_func(t),  m_pose_func(t).rotation_matrix(),
				  m_focal_length, -m_focal_length,
				  -m_detector_origin[0], y - m_detector_origin[1]
				  );
}




class SecondsFrom
{
  boost::posix_time::ptime m_reference;
public:
  inline SecondsFrom( boost::posix_time::ptime const& time ) : m_reference(time) {}

  inline double operator()( boost::posix_time::ptime const& time ) const {
    return double( (time - m_reference).total_microseconds() ) / 1e6;
  }
};

inline boost::posix_time::ptime parse_time(std::string str)
{
  try{
    return boost::posix_time::time_from_string(str);
  }catch(...){
    vw::vw_throw(vw::ArgumentErr() << "Failed to parse time from string: " << str
		 << ". If you are not using Digital Globe images, you may need to specify the session as -t rpc or -t rpcmaprpc.\n");
  }
  return boost::posix_time::time_from_string(str); // Never reached!
}

boost::shared_ptr<DGCameraModel> load_dg_camera_model_from_xml(std::string const& path,
							       bool correct_velocity_aberration)
{
  //vw_out() << "DEBUG - Loading DG camera file: " << camera_file << std::endl;

  // Parse the Digital Globe XML file
  GeometricXML geo;
  AttitudeXML  att;
  EphemerisXML eph;
  ImageXML     img;
  RPCXML       rpc;
  read_xml( path, geo, att, eph, img, rpc );

  // Convert measurements in millimeters to pixels.
  geo.principal_distance /= geo.detector_pixel_pitch;
  geo.detector_origin    /= geo.detector_pixel_pitch;

  // Convert all time measurements to something that boost::date_time can read.
  boost::replace_all( eph.start_time,            "T", " " );
  boost::replace_all( img.tlc_start_time,        "T", " " );
  boost::replace_all( img.first_line_start_time, "T", " " );
  boost::replace_all( att.start_time,            "T", " " );

  // Convert UTC time measurements to line measurements. Ephemeris
  // start time will be our reference frame to calculate seconds against.
  SecondsFrom convert( parse_time( eph.start_time ) );

  // I'm going make the assumption that EPH and ATT are sampled at the same rate and time.
  VW_ASSERT( eph.position_vec.size() == att.quat_vec.size(),
	     vw::MathErr() << "Ephemeris and Attitude don't have the same number of samples." );
  VW_ASSERT( eph.start_time == att.start_time && eph.time_interval == att.time_interval,
	     vw::MathErr() << "Ephemeris and Attitude don't seem to sample with the same t0 or dt." );

  // Convert ephemeris to be position of camera. Change attitude to
  // be the rotation from camera frame to world frame. We also add an
  // additional rotation to the camera frame so X is the horizontal
  // direction to the picture and +Y points down the image (in the direction of flight).
  vw::Quat sensor_coordinate = vw::math::euler_xyz_to_quaternion(vw::Vector3(0,0,geo.detector_rotation - M_PI/2));
  for ( size_t i = 0; i < eph.position_vec.size(); i++ ) {
    eph.position_vec[i] += att.quat_vec[i].rotate( geo.perspective_center );
    att.quat_vec[i] = att.quat_vec[i] * geo.camera_attitude * sensor_coordinate;
  }

  //vw::vw_out() << "DG model load - sensor_coordinate = " << sensor_coordinate << std::endl;
  //geo.printDebugInfo(); // DEBUG INFO

  // Load up the time interpolation class. If the TLCList only has
  // one entry ... then we have to manually drop in the slope and offset.
  if ( img.tlc_vec.size() == 1 ) {
    double direction = 1;
    if ( boost::to_lower_copy( img.scan_direction ) != "forward" ) {
      direction = -1;
    }
    img.tlc_vec.push_back( std::make_pair(img.tlc_vec.front().first +
					  img.avg_line_rate, direction) );
  }

  // Build the TLCTimeInterpolation object and do a quick sanity check.
  vw::camera::TLCTimeInterpolation tlc_time_interpolation( img.tlc_vec,
							       convert( parse_time( img.tlc_start_time ) ) );
  VW_ASSERT( fabs( convert( parse_time( img.first_line_start_time ) ) -
		   tlc_time_interpolation( 0 ) ) < fabs( 1.0 / (10.0 * img.avg_line_rate ) ),
	     vw::MathErr() << "First Line Time and output from TLC lookup table do not agree of the ephemeris time for the first line of the image." );

  vw::Vector2 final_detector_origin
    = subvector(inverse(sensor_coordinate).rotate(vw::Vector3(geo.detector_origin[0],
							      geo.detector_origin[1],
							      0)), 0, 2);

  double et0 = convert( parse_time( eph.start_time ) );
  double at0 = convert( parse_time( att.start_time ) );
  double edt = eph.time_interval;
  double adt = att.time_interval;

  typedef boost::shared_ptr<DGCameraModel> CameraModelPtr;
  return CameraModelPtr(new DGCameraModel(vw::camera::PiecewiseAPositionInterpolation(eph.position_vec, eph.velocity_vec, et0, edt ),
					  vw::camera::LinearPiecewisePositionInterpolation(eph.velocity_vec, et0, edt),
					  vw::camera::SLERPPoseInterpolation(att.quat_vec, at0, adt),
					  tlc_time_interpolation, img.image_size,
					  final_detector_origin,
					  geo.principal_distance, correct_velocity_aberration)
		    );
} // End function load_dg_camera_model()
