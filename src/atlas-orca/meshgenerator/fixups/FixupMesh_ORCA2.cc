
#include "FixupMesh_ORCA2.h"

namespace atlas {
namespace orca {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

void FixupMesh_ORCA2::execute( MeshBuilder& mesh_builder ) const {
    auto mesh = mesh_builder.mesh;

    bool valid_elements{false};
    mesh.metadata().get( "valid_elements", valid_elements );
    if ( valid_elements ) {
        // Nothing to be done here
        return;
    }
    mesh.metadata().set( "valid_elements", true );

    InvalidateElement invalidate( mesh );
    auto& add_element = mesh_builder.add_element;
    auto lonlat       = array::make_view<double, 2>( mesh.nodes().lonlat() );
    auto xy           = array::make_view<double, 2>( mesh.nodes().xy() );

    for ( idx_t n = 0; n < mesh.nodes().size(); ++n ) {
        if ( lonlat( n, LON ) < 70. && lonlat( n, LAT ) < 50. ) {
            lonlat( n, LON ) += 360.;
            xy( n, XX ) = lonlat( n, XX );
        }
    }

    invalidate( 15344 );
    invalidate( 15345 );
    invalidate( 15346 );
    invalidate( 15347 );
    invalidate( 15348 );
    invalidate( 15349 );
    invalidate( 15350 );
    invalidate( 15351 );
    invalidate( 15352 );
    invalidate( 15353 );
    invalidate( 15354 );
    invalidate( 15355 );
    invalidate( 15356 );
    invalidate( 15357 );
    invalidate( 15358 );
    invalidate( 15359 );
    invalidate( 15360 );
    invalidate( 16989 );
    invalidate( 16990 );
    invalidate( 16991 );
    invalidate( 16992 );
    invalidate( 16993 );
    invalidate( 16994 );
    invalidate( 16997 );
    invalidate( 16998 );
    invalidate( 16999 );
    invalidate( 17000 );
    invalidate( 17015 );
    invalidate( 17016 );
    invalidate( 17017 );
    invalidate( 17018 );
    invalidate( 17019 );
    invalidate( 17020 );
    invalidate( 17021 );
    invalidate( 17022 );
    invalidate( 17023 );
    invalidate( 17024 );
    invalidate( 17025 );
    invalidate( 17026 );
    invalidate( 17027 );
    invalidate( 17028 );
    invalidate( 17029 );
    invalidate( 17030 );
    invalidate( 17031 );
    invalidate( 17032 );
    invalidate( 17167 );
    invalidate( 17168 );
    invalidate( 17169 );
    invalidate( 17170 );
    invalidate( 17183 );
    invalidate( 17184 );
    invalidate( 17185 );
    invalidate( 17186 );
    invalidate( 17187 );
    invalidate( 17188 );
    invalidate( 17189 );
    invalidate( 17190 );
    invalidate( 17191 );
    invalidate( 17192 );
    invalidate( 17193 );
    invalidate( 17194 );
    invalidate( 17195 );
    invalidate( 17515 );
    invalidate( 17517 );
    invalidate( 17518 );
    invalidate( 17519 );
    invalidate( 17520 );
    invalidate( 17521 );
    invalidate( 17522 );
    invalidate( 17523 );
    invalidate( 17524 );
    invalidate( 17525 );
    invalidate( 17526 );
    invalidate( 17527 );
    invalidate( 17528 );
    invalidate( 17529 );
    invalidate( 17697 );
    invalidate( 18783 );
    invalidate( 18964 );
    invalidate( 18966 );
    invalidate( 18968 );
    invalidate( 19511 );
    invalidate( 20052 );
    invalidate( 20053 );
    invalidate( 20054 );
    invalidate( 20055 );
    invalidate( 20072 );
    invalidate( 20273 );
    invalidate( 20274 );
    invalidate( 20275 );
    invalidate( 20276 );
    invalidate( 20277 );
    invalidate( 20278 );
    invalidate( 20279 );
    invalidate( 20280 );
    invalidate( 20281 );
    invalidate( 20282 );
    invalidate( 20283 );
    invalidate( 20284 );
    invalidate( 20285 );
    invalidate( 20286 );
    invalidate( 20287 );
    invalidate( 20288 );
    invalidate( 20289 );
    invalidate( 20290 );
    invalidate( 20291 );
    invalidate( 20434 );
    invalidate( 20435 );
    invalidate( 20436 );
    invalidate( 20437 );
    invalidate( 20438 );
    invalidate( 20439 );
    invalidate( 20440 );
    invalidate( 20441 );
    invalidate( 20442 );
    invalidate( 20443 );
    invalidate( 20444 );
    invalidate( 20445 );
    invalidate( 20446 );
    invalidate( 20447 );
    invalidate( 20448 );
    invalidate( 20449 );
    invalidate( 20450 );
    invalidate( 20451 );
    invalidate( 20452 );
    invalidate( 20453 );

    // Other hand-picked invalid elements
    invalidate( 17348 );
    invalidate( 17182 );
    invalidate( 16995 );
    invalidate( 16996 );
    invalidate( 17516 );
    invalidate( 18965 );
    invalidate( 18967 );
    invalidate( 19149 );
    invalidate( 19875 );
    invalidate( 19694 );
    invalidate( 19330 );
    invalidate( 19512 );
    invalidate( 19513 );
    invalidate( 20056 );
    invalidate( 20057 );
    invalidate( 20058 );
    invalidate( 20059 );
    invalidate( 20060 );
    invalidate( 20061 );
    invalidate( 20062 );
    invalidate( 20063 );
    invalidate( 20064 );
    invalidate( 20065 );
    invalidate( 20066 );
    invalidate( 20067 );
    invalidate( 20068 );
    invalidate( 20069 );
    invalidate( 20070 );
    invalidate( 20071 );
    invalidate( 20253 );
    invalidate( 15343 );
    invalidate( 17033 );
    invalidate( 17214 );
    invalidate( 17395 );
    invalidate( 17576 );
    invalidate( 17757 );
    invalidate( 17938 );
    invalidate( 18119 );
    invalidate( 18300 );
    invalidate( 18481 );
    invalidate( 18662 );
    invalidate( 18843 );
    invalidate( 19024 );
    invalidate( 19205 );
    invalidate( 19386 );
    invalidate( 19567 );
    invalidate( 19748 );
    invalidate( 19929 );
    invalidate( 20110 );
    invalidate( 17334 );
    invalidate( 17153 );
    invalidate( 16972 );
    invalidate( 16791 );
    invalidate( 16610 );
    invalidate( 16429 );
    invalidate( 16248 );
    invalidate( 16067 );
    invalidate( 15886 );
    invalidate( 15705 );
    invalidate( 15524 );
    invalidate( 15361 );
    invalidate( 15362 );
    invalidate( 15363 );
    invalidate( 15364 );
    invalidate( 15365 );
    invalidate( 15366 );
    invalidate( 15367 );
    invalidate( 15368 );
    invalidate( 20252 );  // avoid overlap north of caspian sea
    invalidate( 20620 );  // avoid overlap north of caspian sea
    invalidate( 20621 );  // avoid overlap north of caspian sea
    invalidate( 20622 );  // avoid overlap north of caspian sea

    if ( any_of( type, {T, V} ) ) {
        invalidate( 17001 );
    }
    if ( any_of( type, {V, F, U} ) ) {
        invalidate( 18601 );
        invalidate( 18602 );
    }
    if ( any_of( type, {F, U} ) ) {
        invalidate( 17181 );
        invalidate( 20251 );
    }
    if ( any_of( type, {F} ) ) {
        invalidate( 20619 );
    }
    if ( any_of( type, {U} ) ) {
        invalidate( 17696 );
        invalidate( 18782 );
        invalidate( 20251 );
    }

    add_element( 17611, 17794, 17793 );
    if ( any_of( type, {T, F, V} ) ) {
        add_element( 17794, 17795, 17976 );
    }
    if ( any_of( type, {U} ) ) {
        add_element( 17975, 17795, 17976 );
        add_element( 17793, 17794, 17795 );
        add_element( 17793, 17795, 17975 );
    }
    add_element( 17976, 17795, 17977 );
    add_element( 17796, 17795, 17794 );
    add_element( 17797, 17796, 17794 );
    add_element( 17798, 17797, 17794 );
    add_element( 17799, 17798, 17794 );
    add_element( 17800, 17799, 17794 );
    add_element( 17801, 17800, 17794 );
    add_element( 17802, 17801, 17794 );
    add_element( 17803, 17802, 17626 );
    add_element( 17804, 17803, 17626 );
    add_element( 17805, 17804, 17626 );
    add_element( 17806, 17805, 17626 );
    add_element( 17807, 17806, 17626 );
    add_element( 17808, 17807, 17626 );
    add_element( 17626, 17802, 17794 );
    add_element( 17444, 17626, 17794 );
    add_element( 17611, 17444, 17794 );
    add_element( 17447, 17446, 17265 );
    add_element( 17446, 17445, 17265 );
    add_element( 17445, 17444, 17265 );
    add_element( 17445, 17444, 17265 );
    add_element( 17265, 17444, 17611 );
    add_element( 17611, 17429, 17265 );
    add_element( 17429, 17247, 17265 );
    add_element( 17247, 17065, 17265 );
    add_element( 17065, 16883, 17265 );
    add_element( 16883, 16701, 17265 );
    add_element( 16701, 16519, 17265 );
    add_element( 15427, 15428, 15609 );
    add_element( 15428, 15429, 15609 );
    add_element( 15791, 15609, 15429 );
    add_element( 15429, 15430, 15791 );
    add_element( 15430, 15431, 15791 );
    add_element( 15973, 15791, 15431 );
    add_element( 15431, 15432, 15973 );
    add_element( 15432, 15433, 15973 );
    add_element( 16155, 15973, 15433 );
    add_element( 15433, 15434, 16155 );
    add_element( 16337, 16155, 15434 );
    add_element( 15434, 15435, 16337 );
    add_element( 15435, 15436, 16337 );
    add_element( 16519, 16337, 15436 );
    add_element( 15436, 15437, 16519 );
    add_element( 16519, 15437, 15438 );
    add_element( 16519, 15438, 17265 );
    add_element( 15442, 15443, 15610 );
    add_element( 15441, 15442, 15610 );
    add_element( 15610, 15792, 15441 );
    add_element( 15440, 15441, 15792 );
    add_element( 15792, 15974, 15440 );
    add_element( 15439, 15440, 15974 );
    add_element( 15974, 16156, 15439 );
    add_element( 15438, 15439, 16156 );
    add_element( 16156, 17265, 15438 );
    add_element( 17273, 17272, 17066 );
    add_element( 17271, 17270, 17066 );
    add_element( 16884, 17066, 17270 );
    add_element( 17270, 17269, 16884 );
    add_element( 16702, 16884, 17269 );
    add_element( 17269, 17268, 16702 );
    add_element( 16520, 16702, 17268 );
    add_element( 17268, 17267, 16520 );
    add_element( 16338, 16520, 17267 );
    add_element( 17267, 17266, 16338 );
    add_element( 16156, 16338, 17266 );
    add_element( 17266, 17265, 16156 );
    add_element( 15635, 15634, 15453 );
    add_element( 15452, 15453, 15634 );
    add_element( 15634, 15633, 15452 );
    add_element( 15633, 15632, 15452 );
    add_element( 15632, 15631, 15452 );
    add_element( 15451, 15452, 15631 );
    add_element( 15631, 15630, 15451 );
    add_element( 15630, 15629, 15451 );
    add_element( 15629, 15628, 15451 );
    add_element( 15628, 15627, 15451 );
    add_element( 15627, 15626, 15450 );
    add_element( 15626, 15625, 15450 );
    add_element( 15625, 15624, 15450 );
    add_element( 15624, 15623, 15450 );
    add_element( 15623, 15622, 15449 );
    add_element( 15622, 15621, 15449 );
    add_element( 15621, 15620, 15448 );
    add_element( 15620, 15619, 15448 );
    add_element( 15619, 15618, 15447 );
    add_element( 15618, 15617, 15447 );
    add_element( 15617, 15616, 15446 );
    add_element( 15616, 15615, 15446 );
    add_element( 15615, 15614, 15445 );
    add_element( 15614, 15613, 15445 );
    add_element( 15613, 15612, 15444 );
    add_element( 15612, 15611, 15444 );
    add_element( 15611, 15610, 15443 );
    add_element( 15443, 15444, 15611 );
    add_element( 15444, 15445, 15613 );
    add_element( 15445, 15446, 15615 );
    add_element( 15446, 15447, 15617 );
    add_element( 15447, 15448, 15619 );
    add_element( 15448, 15449, 15621 );
    add_element( 15449, 15450, 15623 );
    add_element( 15450, 15451, 15627 );
    add_element( 17277, 17278, 17309 );
    add_element( 17277, 17309, 17308 );
    add_element( 17277, 17308, 17307 );
    if ( any_of( type, {T, V} ) ) {
        add_element( 17094, 17095, 17277 );
        add_element( 17094, 17277, 17307 );
        add_element( 17093, 17094, 17307 );
    }
    if ( any_of( type, {F, U} ) ) {
        add_element( 17093, 17094, 17276 );
        add_element( 17276, 17277, 17307 );
        add_element( 17093, 17276, 17307 );
    }
    add_element( 17093, 17307, 17306 );
    add_element( 17093, 17306, 17305 );
    add_element( 17093, 17305, 17304 );
    add_element( 17093, 17304, 17303 );
    add_element( 17092, 17093, 17303 );
    add_element( 17092, 17303, 17302 );
    add_element( 17092, 17302, 17301 );
    add_element( 17091, 17092, 17301 );
    add_element( 17091, 17301, 17300 );
    add_element( 17090, 17091, 17300 );
    add_element( 17090, 17300, 17299 );
    add_element( 17089, 17090, 17299 );
    add_element( 17088, 17089, 17299 );
    add_element( 17088, 17299, 17298 );
    add_element( 17082, 17083, 17264 );
    add_element( 17261, 17262, 17443 );
    add_element( 17262, 17263, 17625 );
    add_element( 17263, 17264, 17296 );
    add_element( 17262, 17625, 17443 );
    add_element( 17263, 17296, 17625 );
    add_element( 17625, 17296, 17295 );
    add_element( 17625, 17295, 17294 );
    add_element( 17083, 17084, 17264 );
    add_element( 17084, 17085, 17264 );
    add_element( 17264, 17297, 17296 );
    add_element( 17085, 17297, 17264 );
    add_element( 17085, 17086, 17297 );
    add_element( 17086, 17298, 17297 );
    add_element( 17086, 17087, 17298 );
    add_element( 17087, 17088, 17298 );
    add_element( 17624, 17625, 17294 );
    add_element( 17624, 17294, 17293 );
    add_element( 17623, 17624, 17293 );
    add_element( 17622, 17623, 17293 );
    add_element( 17622, 17293, 17292 );
    add_element( 17621, 17622, 17472 );
    add_element( 17472, 17471, 17621 );
    add_element( 17620, 17621, 17471 );
    add_element( 17471, 17470, 17620 );
    add_element( 17619, 17620, 17470 );
    add_element( 17470, 17469, 17619 );
    add_element( 17618, 17619, 17469 );
    add_element( 17469, 17468, 17618 );
    add_element( 17617, 17618, 17468 );
    add_element( 17468, 17467, 17617 );
    add_element( 17616, 17617, 17467 );
    add_element( 17467, 17466, 17616 );
    add_element( 17615, 17616, 17466 );
    add_element( 17466, 17465, 17615 );
    add_element( 17614, 17615, 17465 );
    add_element( 17465, 17464, 17614 );
    add_element( 17613, 17614, 17464 );
    add_element( 17464, 17463, 17613 );
    add_element( 17612, 17613, 17463 );
    add_element( 17463, 17462, 17612 );
    add_element( 17612, 17462, 17461 );
    add_element( 17612, 17461, 17460 );
    if ( any_of( type, {T, V} ) ) {
        add_element( 17276, 17459, 17458 );
        add_element( 17276, 17460, 17459 );
        add_element( 17276, 17612, 17460 );
        add_element( 17430, 17612, 17276 );
        add_element( 17430, 17276, 17275 );
    }
    if ( any_of( type, {F, U} ) ) {
        add_element( 17275, 17459, 17458 );
        add_element( 17275, 17460, 17459 );
        add_element( 17275, 17612, 17460 );
        add_element( 17430, 17612, 17275 );
        add_element( 17275, 17458, 17457 );
    }
    add_element( 17248, 17275, 17274 );
    add_element( 17248, 17274, 17273 );
    add_element( 17066, 17272, 17271 );
    add_element( 17066, 17248, 17273 );
    add_element( 17248, 17430, 17275 );
    if ( any_of( type, {T} ) ) {
        add_element( 18887, 19069, 18886 );
        add_element( 18886, 19069, 19068 );
        add_element( 19069, 19070, 19068 );
    }
    if ( any_of( type, {U} ) ) {
        add_element( 18704, 18705, 18887 );
        add_element( 18703, 18704, 18885 );
        add_element( 18704, 19069, 18885 );
        add_element( 18885, 19069, 19067 );
        add_element( 19069, 19068, 19067 );
        add_element( 19069, 19070, 19068 );
        add_element( 18704, 18887, 19069 );
    }
    if ( any_of( type, {V, F} ) ) {
        add_element( 19069, 19070, 18886 );
        add_element( 18704, 18705, 18887 );
        add_element( 18704, 18887, 19069 );
        add_element( 18703, 18704, 19069 );
        add_element( 18703, 19069, 18885 );
        add_element( 19069, 18886, 18885 );
    }
    if ( any_of( type, {V} ) ) {
        add_element( 19070, 19068, 18886 );
    }
    if ( any_of( type, {F} ) ) {
        add_element( 19070, 19071, 18886 );
        add_element( 18886, 19071, 19068 );
    }
    if ( any_of( type, {T, V, U} ) ) {
        add_element( 19070, 19071, 19068 );
    }
    add_element( 19071, 19437, 19068 );
    add_element( 19071, 19072, 19437 );
    add_element( 19072, 19073, 19255 );
    add_element( 19072, 19255, 19437 );
    add_element( 19068, 19251, 19250 );
    add_element( 19068, 19437, 19251 );
    add_element( 19437, 19619, 19251 );
    add_element( 19619, 19252, 19251 );
    add_element( 19619, 19620, 19252 );
    add_element( 19620, 19253, 19252 );
    add_element( 19620, 19254, 19253 );
    if ( any_of( type, {T, V, F} ) ) {
        add_element( 19620, 19803, 19254 );
        add_element( 19620, 19621, 19803 );
    }
    if ( any_of( type, {U} ) ) {
        add_element( 19620, 19621, 19254 );
        add_element( 19621, 19803, 19254 );
    }
    add_element( 19254, 19803, 19985 );
    add_element( 19254, 19985, 19436 );
    add_element( 19436, 19985, 20167 );
    add_element( 19436, 20167, 19618 );
    add_element( 19618, 20167, 19800 );
    add_element( 20167, 19801, 19800 );
    add_element( 20167, 19802, 19801 );
    add_element( 20168, 19802, 20167 );
    add_element( 20162, 20163, 20344 );
    add_element( 20163, 20164, 20344 );
    add_element( 20164, 20165, 20344 );
    add_element( 20165, 20166, 20344 );
    add_element( 20166, 20345, 20344 );
    add_element( 20166, 19984, 20345 );
    add_element( 20168, 19984, 19802 );
    add_element( 20168, 20345, 19984 );
    add_element( 20168, 20169, 20345 );
    add_element( 20169, 20346, 20345 );
    add_element( 20169, 20170, 20346 );
    add_element( 20170, 20171, 20346 );
    add_element( 20171, 20172, 20346 );
    add_element( 20347, 20346, 20172 );
    add_element( 20172, 20173, 20347 );
    add_element( 20173, 20174, 20347 );
    add_element( 20174, 20348, 20347 );
    add_element( 20174, 20175, 20348 );
    add_element( 20175, 20176, 20348 );
    add_element( 20176, 20349, 20348 );
    add_element( 20176, 20177, 20349 );
    add_element( 20177, 20178, 20349 );
    add_element( 20178, 20179, 20350 );
    add_element( 20178, 20350, 20349 );
    add_element( 20182, 20183, 20365 );
    add_element( 20182, 20365, 20547 );
    add_element( 20179, 20180, 20350 );
    add_element( 20180, 20351, 20350 );
    add_element( 20180, 20181, 20351 );
    add_element( 20181, 20182, 20351 );
    add_element( 20182, 20547, 20351 );
    add_element( 20547, 20352, 20351 );
    add_element( 20547, 20548, 20352 );
    add_element( 20548, 20353, 20352 );
    add_element( 20548, 20549, 20353 );
    add_element( 20549, 20550, 20353 );
    add_element( 20550, 20354, 20353 );
    add_element( 20550, 20551, 20354 );
    add_element( 20551, 20552, 20354 );
    add_element( 20552, 20355, 20354 );
    add_element( 20552, 20553, 20355 );
    add_element( 20553, 20554, 20355 );
    add_element( 20554, 20356, 20355 );
    add_element( 20554, 20555, 20356 );
    add_element( 20555, 20556, 20356 );
    add_element( 20556, 20357, 20356 );
    add_element( 20556, 20557, 20357 );
    add_element( 20557, 20558, 20357 );
    add_element( 20558, 20358, 20357 );
    add_element( 20558, 20559, 20358 );
    add_element( 20559, 20560, 20358 );
    add_element( 20560, 20359, 20358 );
    add_element( 20560, 20561, 20359 );
    add_element( 20561, 20562, 20359 );
    add_element( 20562, 20360, 20359 );
    add_element( 20562, 20563, 20360 );
    add_element( 20563, 20564, 20360 );
    add_element( 20564, 20361, 20360 );
    add_element( 20564, 20565, 20361 );
    add_element( 20565, 20566, 20361 );
    add_element( 20386, 20387, 20362 );
    add_element( 20566, 20362, 20361 );
    add_element( 20387, 20388, 20362 );
    if ( any_of( type, {T, V} ) ) {
        add_element( 20388, 20363, 20362 );
        add_element( 20389, 20390, 20363 );
        add_element( 20390, 20391, 20546 );
        add_element( 20388, 20389, 20363 );
        add_element( 20363, 20390, 20545 );
        add_element( 20545, 20390, 20546 );
    }
    if ( any_of( type, {F, U} ) ) {
        add_element( 20388, 20544, 20362 );
        add_element( 20388, 20389, 20544 );
        add_element( 20389, 20545, 20544 );
        add_element( 20389, 20390, 20545 );
        add_element( 20390, 20391, 20545 );
        add_element( 20391, 20546, 20545 );
    }
    add_element( 20391, 20392, 20546 );
    add_element( 20392, 20728, 20546 );
    add_element( 20392, 20729, 20728 );
    add_element( 20392, 20393, 20729 );
    add_element( 20393, 20730, 20729 );
    add_element( 20393, 20394, 20730 );
    add_element( 20394, 20731, 20730 );
    add_element( 20394, 20395, 20731 );
    add_element( 20395, 20396, 20731 );
    add_element( 20396, 20732, 20731 );
    add_element( 20396, 20397, 20732 );
    add_element( 20397, 20398, 20732 );
    if ( any_of( type, {T, U, V} ) ) {
        add_element( 20398, 20733, 20732 );
        add_element( 20398, 20399, 20733 );
        add_element( 20399, 20400, 20733 );
        add_element( 20400, 20915, 20733 );
    }
    if ( any_of( type, {F} ) ) {
        add_element( 20398, 20914, 20732 );
        add_element( 20398, 20915, 20914 );
        add_element( 20398, 20399, 20915 );
        add_element( 20399, 20400, 20915 );
    }
    add_element( 20400, 20401, 20916 );
    add_element( 20400, 20916, 20915 );
    add_element( 20401, 20402, 20916 );
    add_element( 20402, 20917, 20916 );
    add_element( 20402, 20403, 20917 );
    add_element( 20403, 20918, 20917 );
    add_element( 20403, 20736, 20918 );
    add_element( 20221, 20736, 20403 );
    add_element( 20039, 20736, 20221 );
    add_element( 20039, 20737, 20736 );
    add_element( 20039, 19857, 20737 );
    add_element( 19857, 20738, 20737 );
    add_element( 19857, 19675, 20738 );
    add_element( 17278, 17279, 17309 );
    add_element( 17279, 17280, 17309 );
    add_element( 17491, 17309, 17280 );
    add_element( 17280, 17281, 17491 );
    add_element( 17673, 17491, 17281 );
    add_element( 17855, 17673, 17281 );
    add_element( 18037, 17855, 17281 );
    add_element( 19675, 19493, 20738 );
    add_element( 20739, 20738, 19493 );
    add_element( 19493, 19311, 20739 );
    add_element( 19311, 19129, 20739 );
    add_element( 19129, 18947, 20739 );
    add_element( 18947, 18765, 20739 );
    add_element( 17281, 17282, 18037 );
    add_element( 18219, 18037, 17282 );
    add_element( 18401, 18219, 17282 );
    add_element( 18583, 18401, 17282 );
    add_element( 18765, 18583, 17282 );
    add_element( 20740, 20739, 18765 );
    add_element( 17282, 20740, 18765 );
    add_element( 17282, 20741, 20740 );
    add_element( 17282, 17283, 20741 );
    add_element( 17283, 20742, 20741 );
    add_element( 17283, 17284, 20742 );
    add_element( 17284, 20743, 20742 );
    add_element( 17284, 17285, 20743 );
    add_element( 17285, 20744, 20743 );
    add_element( 17285, 17286, 20744 );
    add_element( 17286, 20745, 20744 );
    add_element( 17286, 17287, 20745 );
    add_element( 17287, 20746, 20745 );
    add_element( 17287, 17288, 20746 );
    add_element( 17288, 20747, 20746 );
    add_element( 17288, 17289, 20747 );
    add_element( 17289, 20748, 20747 );
    add_element( 17289, 17290, 20748 );
    {
        auto& element = add_element( 17109, 20568, 20567 );
        element.ghost = true;  // West ghost element
    }
    {
        auto& element = add_element( 17109, 17110, 20568 );
        element.ghost = true;  // West ghost element
    }
    add_element( 17110, 20569, 20568 );
    add_element( 17110, 17111, 20569 );
    add_element( 17111, 20570, 20569 );
    add_element( 17111, 17112, 20570 );
    add_element( 17112, 20571, 20570 );
    add_element( 17112, 17113, 20571 );
    add_element( 17113, 20572, 20571 );
    add_element( 17113, 17114, 20572 );
    add_element( 17114, 20573, 20572 );
    add_element( 17114, 17115, 20573 );
    add_element( 17115, 20574, 20573 );
    add_element( 17115, 17116, 20574 );
    add_element( 17116, 20575, 20574 );
    add_element( 17116, 17117, 20575 );
    add_element( 17117, 20576, 20575 );
    add_element( 17117, 17118, 20576 );
    add_element( 17118, 20577, 20576 );
    add_element( 17118, 17119, 20577 );
    add_element( 17119, 20578, 20577 );
    add_element( 17119, 17120, 20578 );
    add_element( 17120, 20579, 20578 );
    add_element( 17120, 17121, 20579 );
    add_element( 17121, 20580, 20579 );
    add_element( 17121, 17122, 20580 );
    add_element( 18766, 20580, 17122 );
    add_element( 20404, 20586, 20585 );
    add_element( 20222, 20404, 20585 );
    add_element( 20585, 20584, 20222 );
    add_element( 20040, 20222, 20584 );
    add_element( 20581, 20580, 18766 );
    add_element( 18766, 18948, 20581 );
    add_element( 18948, 19130, 20581 );
    add_element( 20582, 20581, 19130 );
    add_element( 19130, 19312, 20582 );
    add_element( 19312, 19494, 20582 );
    add_element( 20584, 19858, 20040 );
    add_element( 20584, 20583, 19858 );
    add_element( 20583, 20582, 19494 );
    add_element( 19494, 19676, 20583 );
    add_element( 19676, 19858, 20583 );
    add_element( 17122, 17123, 18766 );
    add_element( 17127, 17128, 17310 );
    add_element( 17123, 18584, 18766 );
    add_element( 17123, 18402, 18584 );
    add_element( 17123, 17124, 18402 );
    add_element( 17310, 17492, 17127 );
    add_element( 17126, 17127, 17492 );
    add_element( 17492, 17674, 17126 );
    add_element( 17674, 17856, 17126 );
    add_element( 17125, 17126, 17856 );
    add_element( 17856, 18038, 17125 );
    add_element( 18220, 18402, 17124 );
    add_element( 17125, 18038, 18220 );
    add_element( 17124, 17125, 18220 );
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace orca
}  // namespace atlas
