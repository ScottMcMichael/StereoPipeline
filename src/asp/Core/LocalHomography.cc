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


/// \file LocalHomography.cc
///

#include <vw/Image/ImageView.h>
#include <vw/Image/Transform.h>
#include <vw/Core/ThreadPool.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Stereo/DisparityMap.h>
#include <asp/Core/LocalHomography.h>
#include <asp/Core/StereoSettings.h>
#include <asp/Core/InterestPointMatching.h>

using namespace vw;

namespace asp {

  void split_n_into_k(int n, int k, std::vector<int> & partition){

    VW_ASSERT(n >= k && k > 0,
              ArgumentErr() << "split_n_into_k: Must have n >= k && k > 0.\n");
    int rem = n % k;
    int dx0 = n / k;

    partition.clear();
    int start = 0;
    for (int i = 0; i < k; i++){
      int dx = dx0;
      if (rem > 0){
        dx++;
        rem--;
      }
      partition.push_back(start);
      start += dx;
    }
    partition.push_back(start);

  }

  void write_local_homographies(std::string const& local_hom_file,
                                ImageView<Matrix3x3> const& local_hom){

    std::ofstream fh(local_hom_file.c_str());
    fh.precision(18);
    fh << local_hom.cols() << " " << local_hom.rows() << std::endl;

    for (int col = 0; col < local_hom.cols(); col++){
      for (int row = 0; row < local_hom.rows(); row++){
        Vector<double> V = matrix_to_vector(local_hom(col, row));
        for (int t = 0; t < int(V.size())-1; t++) 
          fh << V[t] << " ";
        if (V.size() > 0) 
          fh << V[V.size()-1] << std::endl;
      }
    }
    fh.close();

    return;
  }

  void read_local_homographies(std::string const& local_hom_file,
                               ImageView<Matrix3x3> & local_hom){

    std::ifstream fh(local_hom_file.c_str());
    if (!fh.good())
      vw_throw( IOErr() << "read_local_homographies: File does not exist: "
                        << local_hom_file << ".\n" );

    int cols, rows;
    if ( !(fh >> cols >> rows) )
      vw_throw( IOErr() << "read_local_homographies: Invalid file: "
                        << local_hom_file << ".\n" );

    local_hom.set_size(cols, rows);
    for (int col = 0; col < local_hom.cols(); col++){
      for (int row = 0; row < local_hom.rows(); row++){

        Vector<double, 9> V;
        for (int t = 0; t < int(V.size()); t++){
          if (! (fh >> V[t]) )
            vw_throw( IOErr() << "read_local_homographies: Invalid file: "
                              << local_hom_file << ".\n" );
        }

        local_hom(col, row) = vector_to_matrix(V);
      }
    }
    fh.close();

    return;
  }

} // namespace asp
