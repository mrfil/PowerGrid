//
//  op_circshift.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/1/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_op_circshift_bones_hpp
#define PowerGrid_op_circshift_bones_hpp


class op_circshift {
public:
    
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_circshift>& in);
    
  template<typename T1> inline static T1 positive_modulo(T1 i, T1 n);
    
};
#endif
