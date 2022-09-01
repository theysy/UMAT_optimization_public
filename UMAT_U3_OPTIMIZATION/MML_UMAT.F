#include "fintrf.h"
      INCLUDE 'LIBRARY/MML_U3SA.FOR'
!     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
!     Declarations
C      implicit none
      IMPLICIT REAL*8(A-H, O-Z)

!     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

!     Function declarations:
      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      integer mxIsNumeric
      mwPointer mxGetM, mxGetN

!     Pointers to input/output mxArrays:
      mwPointer PROPS_PTR, STATEV_PTR, STRESS_PTR, DSTRAN_PTR,
     1          STRESS1_PTR, STATEV1_PTR

!     Array information:
      mwPointer NPROPS, NSTATV, NSTATV2, NTENS, NTENS2
      mwSize PROPS_DIM, STATEV_DIM, TENS_DIM

!     Arguments for computational routine:
      REAL*8, DIMENSION(:), ALLOCATABLE :: PROPS, STATEV, STRESS, 
     1                                       DSTRAN
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: DDSDDE
      COMMON /KUMAT/ NPROPS, NSTATV, NTENS
!-----------------------------------------------------------------------
!     Validate inputs
!     Check that the input is a number.
      if(mxIsNumeric(prhs(1)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:MML_UMAT:NonNumeric',
     1                              'SIG_BAR=NaN')
      endif

!     Get the size of the input array.
      NPROPS = mxGetM(prhs(1))
      NSTATV = mxGetN(prhs(2))
      NTENS = mxGetN(prhs(3))

      ALLOCATE(PROPS(NPROPS), STATEV(NSTATV), STRESS(NTENS),
     1         DSTRAN(NTENS), DDSDDE(NTENS,NTENS))
      PROPS_DIM= NPROPS
      STATEV_DIM= NSTATV
      TENS_DIM= NTENS

!     Create Fortran array from the input argument.
      PROPS_PTR = mxGetPr(prhs(1))
      STATEV_PTR = mxGetPr(prhs(2))
      STRESS_PTR = mxGetPr(prhs(3))
      DSTRAN_PTR = mxGetPr(prhs(4))
      call mxCopyPtrToReal8(PROPS_PTR,PROPS,PROPS_DIM)
      call mxCopyPtrToReal8(STATEV_PTR,STATEV,STATEV_DIM)
      call mxCopyPtrToReal8(STRESS_PTR,STRESS,TENS_DIM)
      call mxCopyPtrToReal8(DSTRAN_PTR,DSTRAN,TENS_DIM)

!     Create matrix for the return argument.
      plhs(1) = mxCreateDoubleMatrix(1,NSTATV,0)
      plhs(2) = mxCreateDoubleMatrix(1,NTENS,0)

      STATEV1_PTR = mxGetPr(plhs(1))
      STRESS1_PTR = mxGetPr(plhs(2))

!     Call the computational subroutine.
      CALL UMAT(PROPS, STATEV, STRESS, DSTRAN, DDSDDE)
      if(ISNAN(STATEV(2))) then
         call mexErrMsgIdAndTxt ('MATLAB:MML_UMAT:NonNumeric',
     1                              'SIG_BAR=NaN')
      endif
!     Load the data into y_ptr, which is the output to MATLAB.
      call mxCopyReal8ToPtr(STRESS,STRESS1_PTR,TENS_DIM)
      call mxCopyReal8ToPtr(STATEV,STATEV1_PTR,STATEV_DIM)
      DEALLOCATE(PROPS, STATEV, STRESS, DSTRAN, DDSDDE)
      return
      end
!-----------------------------------------------------------------------