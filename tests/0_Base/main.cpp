#include "../../operator.h"
int main(int argc, char* argv[]){

  //Assert Correct Index Conversion

  vec m(2);
  m << 3, 3;

  vec dim1(1);  //16 grid
  dim1 << 16;

  vec dim2(2);  //4 x 8 grid
  dim2 << 4, 8;

  vec dim3(3);  //4 x 8 x 16 grid
  dim3 << 4, 8, 16;

  /*
  for(int i = 0; i < dim2.prod(); i++){
    vec n(dim2.size());
    n = op::itop(i, dim2);
    std::cout<<i<<" "<<n<<" "<<op::ptoi(n, dim2)<<std::endl;
    assert(i == op::ptoi(n, dim2));
  }
  */

  for(int i = 0; i < dim3.prod(); i++){
    vec n(dim3.size());
    n = op::itop(i, dim3);
    assert(i == op::ptoi(n, dim3));
  }

  std::cout<<"Index Conversion Assertion Success"<<std::endl;

  //Compute Weights for a 2nd Order Derivate with 3 Support Points!
  std::vector<double> w = op::FD({-1, 0, 1}, 2);

  //1D Matrix Operator
  vec a1(1);
  a1 << -1;
  vec a2(1);
  a2 << 0;
  vec a3(1);
  a3 << 1;

  sparse D1 = op::make({a1, a2, a3}, w, dim1);
  std::cout<<D1<<std::endl;

  //2D Matrix Operator
  vec v1(2);
  v1 << 0, -1;
  vec v2(2);
  v2 << 0, 0;
  vec v3(2);
  v3 << 0, 1;

  sparse D2 = op::make({v1, v2, v3}, w, dim2);
  std::cout<<D2<<std::endl;

  return 0;
}
