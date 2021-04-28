/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>


#include "atlas-orca/util/Enums.h"
#include "atlas-orca/util/Flag.h"
#include "atlas-orca/util/OrcaData.h"
#include "atlas-orca/util/OrcaDataFile.h"
#include "atlas-orca/util/OrcaPeriodicity.h"
#include "atlas/runtime/Exception.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/log/ProgressTimer.h"

namespace atlas {
namespace orca {


struct ReadLine {
    bool verbose = false;
    int version  = 2;
    double lat;
    double lon;
    double water;
    double core;
    int I;
    int J;
    int I_master;
    int J_master;
    double lat_master;
    double lon_master;
    friend std::istream& operator>>( std::istream& in, ReadLine& r ) {
        if ( r.version >= 1 ) {
            in >> r.lat;
            in >> r.lon;
            in >> r.water;
            in >> r.core;
        }
        if ( r.version >= 2 ) {
            in >> r.I;
            in >> r.J;
            in >> r.I_master;
            in >> r.J_master;
            in >> r.lat_master;
            in >> r.lon_master;
        }
        return in;
    }

    void fix() {
        if ( version == 1 ) {
            return;
        }
        if ( core == 0. ) {
            if ( I == I_master && J == J_master ) {
                if ( verbose ) {
                    Log::info() << "Point {i,j} = {" << I - 1 << "," << J - 1 << "} --> {" << I_master - 1 << ","
                                << J_master - 1 << "} should be marked as core" << std::endl;
                }
                core = 1.;
            }
        }
        else {
            if ( I != I_master || J != J_master ) {
                if ( verbose ) {
                    Log::info() << "Point {I,J} = {" << I << "," << J << "} --> {" << I_master << "," << J_master
                                << "} should be marked as non-core" << std::endl;
                }
                core = 0.;
            }
        }
        if ( lat != lat_master && J_master != 0 ) {
            if ( verbose ) {
                Log::info() << "Point {I,J} = {" << I << "," << J << "} --> {" << I_master << "," << J_master
                            << "} has wrong lat: [" << lat << " != " << lat_master << "]" << std::endl;
            }
            lat = lat_master;
        }
        if ( lon != lon_master && J_master != 0 ) {
            if ( verbose ) {
                Log::info() << "Point {I,J} = {" << I << "," << J << "} --> {" << I_master << "," << J_master
                            << "} has wrong lon: [" << lon << " != " << lon_master << "]" << std::endl;
            }
            lon = lon_master;
        }
    }
};

class AsciiReader {
private:
    bool verbose_{false};
    int version_{2};
    std::vector<double> pivot_v1_;

public:
    AsciiReader( const util::Config& config ) {
        config.get( "verbose", verbose_ );
        config.get( "version", version_ );
        config.get( "pivot", pivot_v1_ );
    }

    void read( const std::string& uri, OrcaData& data ) {
        OrcaDataFile file{uri};

        auto trace = atlas::Trace( Here(), "read" );

        std::ifstream ifstrm{file.c_str()};

        // Reading header
        std::string line;
        std::getline( ifstrm, line );
        std::istringstream iss{line};

        std::int32_t nx_halo = 0;
        std::int32_t ny_halo = 0;
        auto& halo           = data.halo;
        halo                 = {0, 0, 0, 0};

        ATLAS_ASSERT( iss >> nx_halo >> ny_halo, "Error while reading header" );
        if ( version_ == 1 ) {
            int periodicity;
            iss >> periodicity >> halo[HALO_EAST] >> halo[HALO_WEST] >> halo[HALO_SOUTH] >> halo[HALO_NORTH];
        }

        size_t jstride = nx_halo;

        size_t size = nx_halo * ny_halo;

        data.lon.resize( size );
        data.lat.resize( size );
        data.flags.resize( size );


        std::vector<idx_t> master( size );

        // Reading coordinates

        {
            eckit::Channel blackhole;
            eckit::ProgressTimer progress( "Reading file " + file.path(), size, " point", double( 1 ),
                                           size > 5.e6 && verbose_ ? Log::info() : blackhole );


            ReadLine r;
            r.version = version_;
            r.verbose = verbose_;


            for ( size_t n = 0; n < size; ++n ) {
                ++progress;
                master[n] = n;
                std::getline( ifstrm, line );
                std::istringstream iss{line};
                iss >> r;
                r.fix();

                data.lon[n] = r.lon;
                data.lat[n] = r.lat;

                Flag flag{data.flags[n]};
                if ( r.core == 0. ) {  // ghost
                    flag.set( Flag::GHOST );
                    if ( version_ >= 2 ) {
                        master[n] = ( r.I_master - 1 ) + ( r.J_master - 1 ) * std::int64_t( jstride );
                        if ( master[n] < 0 ) {
                            master[n] = -1;
                        }
                    }
                }
                if ( r.water != 0. ) {
                    flag.set( Flag::WATER );
                }
                if ( r.J == ny_halo && r.I == nx_halo / 2 + 2 && version_ >= 2 ) {
                    ATLAS_ASSERT( r.I != r.I_master );
                    ATLAS_ASSERT( r.J != r.J_master );
                    data.pivot[0] = double( ( r.I - 1 ) + ( r.I_master - 1 ) ) * 0.5;
                    data.pivot[1] = double( ( r.J - 1 ) + ( r.J_master - 1 ) ) * 0.5;
                }
            }
            if ( version_ < 2 ) {
                ATLAS_ASSERT( pivot_v1_.size() == 2 );
                std::copy( pivot_v1_.begin(), pivot_v1_.end(), data.pivot.begin() );
            }
            ATLAS_ASSERT( data.pivot[0] >= 0. );
            ATLAS_ASSERT( data.pivot[1] >= 0. );
        }

        trace.stop();
        if ( trace.elapsed() > 1. && verbose_ ) {
            Log::info() << "Reading file " << file.path() << " took " << trace.elapsed() << " seconds." << std::endl;
        }

        if ( version_ >= 2 ) {
            auto is_ghost = [&]( idx_t i, idx_t j ) {
                gidx_t n = j * nx_halo + i;
                return n != master[n];
            };
            auto j_contains_core_point = [&]( idx_t j ) {
                bool contains_core_point = false;
                for ( idx_t i = 0; i < nx_halo; ++i ) {
                    if ( not is_ghost( i, j ) ) {
                        contains_core_point = true;
                        break;
                    }
                }
                return contains_core_point;
            };
            auto i_contains_core_point = [&]( idx_t i ) {
                bool contains_core_point = false;
                for ( idx_t j = halo[HALO_SOUTH]; j < ny_halo - halo[HALO_NORTH]; ++j ) {
                    if ( not is_ghost( i, j ) ) {
                        contains_core_point = true;
                        break;
                    }
                }
                return contains_core_point;
            };

            for ( idx_t j = 0; j < ny_halo; ++j, ++halo[HALO_SOUTH] ) {
                if ( j_contains_core_point( j ) ) {
                    break;
                }
            }
            for ( idx_t j = ny_halo - 1; j >= 0; --j, ++halo[HALO_NORTH] ) {
                if ( j_contains_core_point( j ) ) {
                    break;
                }
            }
            for ( idx_t i = 0; i < nx_halo; ++i, ++halo[HALO_WEST] ) {
                if ( i_contains_core_point( i ) ) {
                    break;
                }
            }
            for ( idx_t i = nx_halo - 1; i >= 0; --i, ++halo[HALO_EAST] ) {
                if ( i_contains_core_point( i ) ) {
                    break;
                }
            }
        }
        data.dimensions[0] = nx_halo;
        data.dimensions[1] = ny_halo;

        if ( version_ >= 2 ) {
            validate( data, master );
        }
        data.setGhost();
        data.makeHaloConsistent();
    }


    void validate( const OrcaData& data, const std::vector<idx_t>& check_master ) {
        constexpr int N = HALO_NORTH;
        constexpr int W = HALO_WEST;
        constexpr int S = HALO_SOUTH;
        constexpr int E = HALO_EAST;
        auto& halo      = data.halo;

        OrcaPeriodicity compute_master{data};

        Log::info() << "pivot_i = " << data.pivot[0] << std::endl;
        Log::info() << "pivot_j = " << data.pivot[1] << std::endl;
        Log::info() << "---------------------------------------------------" << std::endl;
        idx_t ni = data.dimensions[0];
        idx_t nj = data.dimensions[1];
        for ( idx_t j = 0; j < nj; ++j ) {
            for ( idx_t i = 0; i < ni; ++i ) {
                auto master = compute_master( i, j );
                idx_t n     = ni * j + i;
                Flag flag{data.flags[n]};

                if ( j < halo[S] || j >= nj - halo[N] ) {
                    ATLAS_ASSERT( flag.test( Flag::GHOST ) );
                }
                idx_t n_master = master.j * ni + master.i;
                if ( j > halo[S] ) {
                    if ( flag.test( Flag::GHOST ) != ( n_master != n ) ) {
                        Log::warning() << "Mismatch in expected ghost:" << std::endl;
                        Log::warning() << "    point  " << n << "\t" << PointIJ( i, j ) << std::endl;
                        Log::warning() << "    master " << n_master << "\t" << master << std::endl;
                        Log::warning() << "    ghost: " << std::boolalpha << flag.test( Flag::GHOST ) << std::endl;
                        Log::warning() << "---------------------------------------------------" << std::endl;
                    }

                    if ( n_master != check_master[n] ) {
                        idx_t i_check = check_master[n] % ni;
                        idx_t j_check = check_master[n] / ni;
                        if ( check_master[n] < 0 ) {
                            Log::warning() << "Ignored (1) : ";  //continue;
                        }
                        else if ( ( double( j_check ) > data.pivot[1] ) ) {
                            Log::warning() << "Ignored (2) : ";  //continue;
                        }
                        else if ( ( double( j_check ) == data.pivot[1] ) && ( i_check >= ni - halo[E] - halo[W] ) ) {
                            Log::warning() << "Ignored (3) : ";  //continue;
                        }
                        Log::warning() << "Mismatch in master index:" << std::endl;
                        Log::warning() << "    point  " << n << "\t" << PointIJ( i, j ) << std::endl;
                        Log::warning() << "    master " << n_master << "\t" << master << std::endl;
                        Log::warning() << "    file   " << check_master[n] << "\t" << PointIJ( i_check, j_check )
                                       << std::endl;
                        Log::warning() << "---------------------------------------------------" << std::endl;
                        //ATLAS_ASSERT( n_master == master[n] );
                    }
                }
            }
        }
    }
};

}  // namespace orca
}  // namespace atlas

//------------------------------------------------------------------------------------------------------
