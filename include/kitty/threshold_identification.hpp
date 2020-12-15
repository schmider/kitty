/* kitty: C++ truth table library
 * Copyright (C) 2017-2020  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file threshold_identification.hpp
  \brief Threshold logic function identification

  \author CS-472 2020 Fall students
*/

#pragma once

#include <vector>
#include <lpsolve/lp_lib.h> /* uncomment this line to include lp_solve */
#include "traits.hpp"
#include "print.hpp"
#include "cube.hpp"
#include "bit_operations.hpp"
#include "properties.hpp"
#include "isop.hpp"
#include "implicant.hpp"
#include "operations.hpp"

namespace kitty
{

/*! \brief Threshold logic function identification

  Given a truth table, this function determines whether it is a threshold logic function (TF)
  and finds a linear form if it is. A Boolean function is a TF if it can be expressed as

  f(x_1, ..., x_n) = \sum_{i=1}^n w_i x_i >= T

  where w_i are the weight values and T is the threshold value.
  The linear form of a TF is the vector [w_1, ..., w_n; T].

  \param tt The truth table
  \param plf Pointer to a vector that will hold a linear form of `tt` if it is a TF.
             The linear form has `tt.num_vars()` weight values and the threshold value
             in the end.
  \return `true` if `tt` is a TF; `false` if `tt` is a non-TF.
*/

template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
inline bool less( const TT& first, const TT& second )
{
  if ( first.num_vars() != second.num_vars() )
  {
    return false;
  }
    for(int i = 0;  i <first.num_bits(); i++){
        if (!(get_bit(first,i) <= get_bit(second,i))) {
            return false;
        }
    }
    return true;
}
template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
inline bool greater( const TT& first, const TT& second )
{
    if ( first.num_vars() != second.num_vars() )
    {
      return false;
    }
    for(int i = 0;  i <first.num_bits(); i++){
        if (!(get_bit(first,i) >= get_bit(second,i))) {
            return false;
        }
    }
    return true;
  }
template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
bool is_threshold( const TT& tt, std::vector<int64_t>* plf = nullptr )
{
    std::vector<int64_t> modified_variable (tt.num_vars(),0);
  std::vector<int64_t> linear_form;
    auto positive_unate = tt;
    std::vector<int64_t> positive_unate_sol;
  /* TODO */
  /* if tt is non-TF: */
    
    for(int i = 0; i < positive_unate.num_vars(); i++){
       /* std::cout <<" tonio cofactor1 " << to_binary(cofactor1(positive_unate, i)) << "\n";
        std::cout <<" tonio cofactor0 " << to_binary(cofactor0(positive_unate, i)) << "\n";
        std::cout <<" tonio less " << less(cofactor1(positive_unate, i),cofactor0(positive_unate, i)) << "\n";
        std::cout <<" tonio more " << greater(cofactor1(positive_unate, i),cofactor0(positive_unate, i)) << "\n";*/
        if(less(cofactor1(positive_unate, i),cofactor0(positive_unate, i))){
            kitty::flip_inplace(positive_unate, i);
            modified_variable.at(i) = 1;
        }else if(greater(cofactor1(positive_unate, i),cofactor0(positive_unate, i))){
            
        }else{
            return false;
        }
    }
   const auto cubes =  isop(positive_unate);
 
     for ( auto cube : cubes )
     {
       cube.print( tt.num_vars() );
       std::cout << std::endl;
         
     }
    const auto cubes_neg =  isop(unary_not(positive_unate));
  
      for ( auto cube : cubes_neg )
      {
        cube.print( tt.num_vars() );
        std::cout << std::endl;
          
      }
    
      lprec *lp;
      int Ncol, *colno = NULL, ret = 0;
      REAL *row = NULL;

      /* We will build the model row by row
         So we start with creating a model with 0 rows and 2 columns */
      Ncol = tt.num_vars() + 1;
      lp = make_lp(0, Ncol);
      if(lp == NULL)
        ret = 1; /* couldn't construct a new model... */

      if(ret == 0) {
        /* let us name our variables. Not required, but can be useful for debugging */
          const char* letter[26] = {
              "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M",
              "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"};
          for(int i = 1; i <= tt.num_vars(); i++){
              set_col_name(lp, i, (char*)letter[i-1]);
          }
          set_col_name(lp, tt.num_vars() + 1, "T");
          
        /* create space large enough for one row */
        colno = (int *) malloc(Ncol * sizeof(*colno));
        row = (REAL *) malloc(Ncol * sizeof(*row));
        if((colno == NULL) || (row == NULL))
          ret = 2;
      }

    if(ret == 0) {
      set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */
        
        for(int j = 0; j <= tt.num_vars(); j++){
            colno[0] = j+1;
            row[j] = 1;
            add_constraintex(lp, 1, row, colno, GE, 0);
        }
        
    
        for(int i = 0; i < cubes.size() ; i++){
            auto cube = cubes.at(i);
            for(int j = 0; j < tt.num_vars(); j++){
                if(cube.get_mask(j)==0){
                    colno[j] = j+1;
                    row[j] = 0;
                }else{
                    if(cube.get_bit(j)==1){
                        colno[j] = j+1;
                        row[j] = 1;
                    }else{
                        colno[j] = j+1;
                        row[j] = 0;
                    }
                }
            }
            colno[tt.num_vars()] = tt.num_vars() + 1;
            row[tt.num_vars()] = -1;
            add_constraintex(lp, tt.num_vars()+1, row, colno, GE, 0);
        }
      }
    if(ret == 0) {
      
        for(int i = 0; i < cubes_neg.size() ; i++){
            auto cube = cubes_neg.at(i);
            cube.print( tt.num_vars() );
            std::cout << std::endl;
            for(int j = 0; j < tt.num_vars(); j++){
                if(cube.get_mask(j)==0){
                    colno[j] = j+1;
                    row[j] = 1;
                }else{
                    if(cube.get_bit(j)==0){
                        colno[j] = j+1;
                        row[j] = 0;
                    }else{
                        colno[j] = j+1;
                        row[j] = 1;
                    }
                }
                
            }
            colno[tt.num_vars()] = tt.num_vars() + 1;
            row[tt.num_vars()] = -1;
            add_constraintex(lp, tt.num_vars() + 1, row, colno, LE, -1);
        }
      }


      if(ret == 0) {
        set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */
          int j = 0;
          for(int i = 1; i <= tt.num_vars(); i++){
                  colno[j] = i;
                  row[j++] = 1;
          }
          colno[j] = tt.num_vars() + 1;
          row[j++] = 1;
          set_obj_fnex(lp, j, row, colno);
      }

      if(ret == 0) {
        /* set the object direction to maximize */
        set_minim(lp);

        /* just out of curioucity, now show the model in lp format on screen */
        /* this only works if this is a console application. If not, use write_lp and a filename */
        write_LP(lp, stdout);
        /* write_lp(lp, "model.lp"); */


        /* Now let lpsolve calculate a solution */
        ret = solve(lp);
          if(ret != 0){
              return false;
              
          }
      }

      if(ret == 0) {
        /* a solution is calculated, now lets get some results */

        /* objective value */
        printf("Objective value: %f\n", get_objective(lp));

        /* variable values */
        get_variables(lp, row);
        for(int j = 0; j < Ncol; j++){
          printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);
            positive_unate_sol.push_back(row[j]);
         }
        /* we are done now */
      }
    for(int i = 0; i< tt.num_vars(); i++){
        if(modified_variable.at(i)==1){
            int dummy = positive_unate_sol.at(i);
            positive_unate_sol.at(i) = -dummy;
            positive_unate_sol.at(tt.num_vars()) += -dummy;
        }
    }
    linear_form = positive_unate_sol;
      /* free allocated memory */
      if(row != NULL)
        free(row);
      if(colno != NULL)
        free(colno);

      if(lp != NULL) {
        /* clean up such that all used memory by lpsolve is freed */
        delete_lp(lp);
      }


  /* if tt is TF: */
  /* push the weight and threshold values into `linear_form` */
  if ( plf )
  {
    *plf = linear_form;
  }
  return true;
}

} /* namespace kitty */
