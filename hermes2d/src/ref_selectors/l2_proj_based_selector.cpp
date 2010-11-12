#include "../h2d_common.h"
#include "../../../hermes_common/matrix.h"
#include "../solution.h"
#include "../shapeset/shapeset_l2_all.h"
#include "../element_to_refine.h"
#include "l2_proj_based_selector.h"

namespace RefinementSelectors {
  L2Shapeset L2ProjBasedSelector::default_shapeset;

  const int L2ProjBasedSelector::H2DRS_MAX_L2_ORDER = H2DRS_MAX_ORDER;

  L2ProjBasedSelector::L2ProjBasedSelector(CandList cand_list, double conv_exp, int max_order, L2Shapeset* user_shapeset)
    : ProjBasedSelector(cand_list, conv_exp, max_order, user_shapeset == NULL ? &default_shapeset : user_shapeset, Range<int>(1,1), Range<int>(0, H2DRS_MAX_L2_ORDER)) {}

  void L2ProjBasedSelector::set_current_order_range(Element* element) {
    current_max_order = this->max_order;
    if (current_max_order == H2DRS_DEFAULT_ORDER)
      current_max_order = (20 - element->iro_cache)/2 - 2; // default
    else
      current_max_order = std::min(current_max_order, (20 - element->iro_cache)/2 - 2); // user specified
    current_min_order = 1;
  }

  void L2ProjBasedSelector::precalc_shapes(const double3* gip_points, const int num_gip_points, const Trf* trfs, const int num_noni_trfs, const std::vector<ShapeInx>& shapes, const int max_shape_inx, TrfShape& svals) {
    //for all transformations
    bool done = false;
    int inx_trf = 0;
    while (!done && inx_trf < H2D_TRF_NUM) {
      //prepare data for processing
      const Trf& trf = trfs[inx_trf];
      std::vector<TrfShapeExp>& trf_svals = svals[inx_trf];

      //allocate
      trf_svals.resize(max_shape_inx + 1);

      //for all shapes
      const int num_shapes = (int)shapes.size();
      for(int i = 0; i < num_shapes; i++) {
        int inx_shape = shapes[i].inx;
        TrfShapeExp& shape_exp = trf_svals[inx_shape];

        //allocate
        shape_exp.allocate(H2D_L2FE_NUM, num_gip_points);

        //for all GIP points
        for(int k = 0; k < num_gip_points; k++) {
          //transform coordinates
          double ref_x = gip_points[k][H2D_GIP2D_X] * trf.m[0] + trf.t[0];
          double ref_y = gip_points[k][H2D_GIP2D_Y] * trf.m[1] + trf.t[1];

          //for all expansions: retrieve values
          shape_exp[H2D_L2FE_VALUE][k] = shapeset->get_fn_value(inx_shape, ref_x, ref_y, 0);
        }
      }

      //move to the next transformation
      if (inx_trf == H2D_TRF_IDENTITY)
        done = true;
      else {
        inx_trf++;
        if (inx_trf >= num_noni_trfs) //if all transformations were processed, move to the identity transformation
          inx_trf = H2D_TRF_IDENTITY;
      }
    }
    error_if(!done, "All transformation processed but identity transformation not found."); //identity transformation has to be the last transformation
  }

  void L2ProjBasedSelector::precalc_ortho_shapes(const double3* gip_points, const int num_gip_points, const Trf* trfs, const int num_noni_trfs, const std::vector<ShapeInx>& shapes, const int max_shape_inx, TrfShape& svals) {
    //calculate values
    precalc_shapes(gip_points, num_gip_points, trfs, num_noni_trfs, shapes, max_shape_inx, svals);

    //calculate orthonormal basis
    const int num_shapes = (int)shapes.size();
    for(int i = 0; i < num_shapes; i++) {
      const int inx_shape_i = shapes[i].inx;

      //orthogonalize
      for(int j = 0; j < i; j++) {
        const int inx_shape_j = shapes[j].inx;

        //calculate product of non-transformed functions
        double product = 0.0;
        for(int k = 0; k < num_gip_points; k++) {
          double sum = 0.0;
          sum += svals[H2D_TRF_IDENTITY][inx_shape_i][H2D_L2FE_VALUE][k] * svals[H2D_TRF_IDENTITY][inx_shape_j][H2D_L2FE_VALUE][k];
          product += gip_points[k][H2D_GIP2D_W] * sum;
        }

        //for all transformations
        int inx_trf = 0;
        bool done = false;
        while (!done && inx_trf < H2D_TRF_NUM) {
          //for all integration points
          for(int k = 0; k < num_gip_points; k++) {
            svals[inx_trf][inx_shape_i][H2D_L2FE_VALUE][k] -= product * svals[inx_trf][inx_shape_j][H2D_L2FE_VALUE][k];
          }

          //move to the next transformation
          if (inx_trf == H2D_TRF_IDENTITY)
            done = true;
          else {
            inx_trf++;
            if (inx_trf >= num_noni_trfs) //if all transformations were processed, move to the identity transformation
              inx_trf = H2D_TRF_IDENTITY;
          }
        }
        error_if(!done, "All transformation processed but identity transformation not found."); //identity transformation has to be the last transformation
      }

      //normalize
      //calculate norm
      double norm_squared = 0.0;
      for(int k = 0; k < num_gip_points; k++) {
        double sum = 0.0;
        sum += sqr(svals[H2D_TRF_IDENTITY][inx_shape_i][H2D_L2FE_VALUE][k]);
        norm_squared += gip_points[k][H2D_GIP2D_W] * sum;
      }
      double norm = sqrt(norm_squared);
      assert_msg(finite(1/norm), "Norm (%g) is almost zero.", norm);

      //for all transformations: normalize
      int inx_trf = 0;
      bool done = false;
      while (!done && inx_trf < H2D_TRF_NUM) {
        //for all integration points
        for(int k = 0; k < num_gip_points; k++) {
          svals[inx_trf][inx_shape_i][H2D_L2FE_VALUE][k] /= norm;
        }

        //move to the next transformation
        if (inx_trf == H2D_TRF_IDENTITY)
          done = true;
        else {
          inx_trf++;
          if (inx_trf >= num_noni_trfs) //if all transformations were processed, move to the identity transformation
            inx_trf = H2D_TRF_IDENTITY;
        }
      }
      error_if(!done, "All transformation processed but identity transformation not found."); //identity transformation has to be the last transformation
    }
  }

  scalar** L2ProjBasedSelector::precalc_ref_solution(int inx_son, Solution* rsln, Element* element, int intr_gip_order) {
    //set element and integration order
    rsln->set_active_element(element);
    rsln->set_quad_order(intr_gip_order);

    //fill with values
    scalar** rvals_son = precalc_rvals[inx_son];
    rvals_son[H2D_L2FE_VALUE] = rsln->get_fn_values(0);

    return rvals_son;
  }

  double** L2ProjBasedSelector::build_projection_matrix(double3* gip_points, int num_gip_points,
    const int* shape_inx, const int num_shapes) {
    //allocate
    double** matrix = new_matrix<double>(num_shapes, num_shapes);

    //calculate products
    int inx_row = 0;
    for(int i = 0; i < num_shapes; i++, inx_row += num_shapes) {
      double* matrix_row = matrix[i];
      int shape0_inx = shape_inx[i];
      for(int k = 0; k < num_shapes; k++) {
        int shape1_inx = shape_inx[k];

        double value = 0.0;
        for(int j = 0; j < num_gip_points; j++) {
          double gip_x = gip_points[j][H2D_GIP2D_X], gip_y = gip_points[j][H2D_GIP2D_Y];
          double value0 = shapeset->get_value(H2D_FEI_VALUE, shape0_inx, gip_x, gip_y, 0);
          double value1 = shapeset->get_value(H2D_FEI_VALUE, shape1_inx, gip_x, gip_y, 0);

          value += gip_points[j][H2D_GIP2D_W] * (value0*value1);
        }

        matrix_row[k] = value;
      }
    }

    return matrix;
  }

  scalar L2ProjBasedSelector::evaluate_rhs_subdomain(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, const ElemSubShapeFunc& sub_shape) {
    scalar total_value = 0;
    for(int gip_inx = 0; gip_inx < sub_gip.num_gip_points; gip_inx++) {
      //get location and transform it
      double3 &gip_pt = sub_gip.gip_points[gip_inx];

      //get value of a shape function
      double shape_value = sub_shape.svals[H2D_L2FE_VALUE][gip_inx];

      ////DEBUG-BEGIN
      //double ref_x = gip_pt[H2D_GIP2D_X] * sub_trf.trf->m[0] + sub_trf.trf->t[0];
      //double ref_y = gip_pt[H2D_GIP2D_Y] * sub_trf.trf->m[1] + sub_trf.trf->t[1];
      //double shape_valueA = shapeset->get_fn_value(sub_shape.inx, ref_x, ref_y, 0);
      //error_if(std::abs(shape_value - shape_valueA) > 1E-14, "A0");
      ////DEBUG-END

      //get value of ref. solution
      scalar ref_value;
      ref_value = sub_gip.rvals[H2D_L2FE_VALUE][gip_inx];

      //evaluate a right-hand value
      scalar value = (shape_value * ref_value);

      total_value += gip_pt[H2D_GIP2D_W] * value;
    }
    return total_value;
  }

  double L2ProjBasedSelector::evaluate_error_squared_subdomain(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, const ElemProj& elem_proj) {
    double total_error_squared = 0;
    for(int gip_inx = 0; gip_inx < sub_gip.num_gip_points; gip_inx++) {
      //get location and transform it
      double3 &gip_pt = sub_gip.gip_points[gip_inx];

      //calculate value of projected solution
      scalar proj_value = 0;
      for(int i = 0; i < elem_proj.num_shapes; i++) {
        int shape_inx = elem_proj.shape_inxs[i];
        proj_value += elem_proj.shape_coefs[i] * elem_proj.svals[shape_inx][H2D_L2FE_VALUE][gip_inx];
      }

      //get value of ref. solution
      scalar ref_value = sub_gip.rvals[H2D_L2FE_VALUE][gip_inx];

      //evaluate error
      double error_squared = sqr(proj_value - ref_value);

      total_error_squared += gip_pt[H2D_GIP2D_W] * error_squared;
    }
    return total_error_squared;
  }
}

