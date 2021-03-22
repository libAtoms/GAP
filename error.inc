! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   GAP (Gaussian Approximation Potental)
! HND X   
! HND X
! HND X   Portions of GAP were written by Albert Bartok-Partay, Gabor Csanyi, 
! HND X   Copyright 2006-2021.
! HND X
! HND X   Portions of GAP were written by Noam Bernstein as part of
! HND X   his employment for the U.S. Government, and are not subject
! HND X   to copyright in the USA.
! HND X
! HND X   GAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   GAP is distributed in the hope that it will be useful for non-commercial
! HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied 
! HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! HND X   ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original licensors,
! HND X   Gabor Csanyi or Albert Bartok-Partay. The ASL is also published at
! HND X   http://github.com/gabor1/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   A. P. Bartok et al Physical Review Letters vol 104 p136403 (2010)
! HND X
! HND X   When using the SOAP kernel or its variants, please additionally cite:
! HND X
! HND X   A. P. Bartok et al Physical Review B vol 87 p184115 (2013)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
!X 
!X Error handling, see error.f95 for the functions called in these macros. 
!X 
!X Error passing works as follows: 
!X - *error* needs to be intent(out) and optional 
!X - all functions that receive *error* as an argument must call INIT_ERROR(error) 
!X - RAISE_ERROR is used whenever an error occurs. If *error* is not present, 
!X   the program execution will be terminated immediately. If *error* is 
!X   present it will be set to some value not equal ERROR_NONE and the execution 
!X   of the subroutine will be stopped. 
!X - PASS_ERROR is used after a function or subroutine that returns error, i.e. 
!X     call sub(..., error=error) 
!X     PASS_ERROR(error) 
!X   If no error occurs (i.e. error==ERROR_NONE), execution will proceed as 
!X   usual. If an error occured, the current function will be terminated after 
!X   the location of the error is passed to the error module. 
!X   If the calling routine handles the error itself, rather than passing 
!X   it up with PASS_ERROR(), CLEAR_ERROR() should be used to clear the error 
!X   info stack 
!X - PASS_ERROR_WITH_INFO is like PASS_ERROR, just an additional string can be 
!X   provided describing the error, or parameters. 
!X - HANDLE_ERROR will print the error history and stop execution of the program 
!X   after an error occured. 
!X 
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
  
#define INIT_ERROR(error)  if (present(error)) then ; error = ERROR_NONE ; endif 
#define ASSERT(condition, message, error)   if (.not. (condition)) then ; RAISE_ERROR(message, error) ; endif 
  
#define RAISE_ERROR(message, error)   if (.true.) then ; call push_error_with_info(message, __FILE__, __LINE__) ; if (present(error)) then ; error = ERROR_UNSPECIFIED ; return ; else ; call error_abort(error) ; endif ; endif 
  
#define RAISE_ERROR_WITH_KIND(kind, message, error)   if (.true.) then ; call push_error_with_info(message, __FILE__, __LINE__, kind) ; if (present(error)) then ; error = kind ; return ; else ; call error_abort(error) ; endif ; endif 
  
#define PASS_ERROR(error)  if (present(error)) then ; if (error /= ERROR_NONE) then ; call push_error(__FILE__, __LINE__) ; return ; endif ; endif 
  
#define PASS_ERROR_WITH_INFO(message, error)  if (present(error)) then ; if (error /= ERROR_NONE) then ; call push_error_with_info(message, __FILE__, __LINE__) ; return ; endif ; endif 
  
#define HANDLE_ERROR(error)  if (error /= ERROR_NONE) then ; call push_error(__FILE__, __LINE__) ; call error_abort(error) ; endif 
  
#define CLEAR_ERROR(error) call error_clear_stack() 
   
#define PRINT_LINE_NUMBER if(.true.) then; print "('LINE ' i0)",__LINE__; endif  
  
  
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
!X  
!X MPI errors  
!X  
!X MPI error string are obtained using mpi_error_string and then pushed  
!X onto the error stack.  
!X  
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
  
#define PASS_MPI_ERROR(mperror, error)  if (mperror /= MPI_SUCCESS) then ; call push_MPI_error(mperror, __FILE__, __LINE__) ; if (present(error)) then ; error = ERROR_MPI ; return ; else ; call error_abort(error) ; endif ; endif  
  
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
!X  
!X MPI BCAST errors  
!X  
!X Extension of error handling macros to MPI cases where processes  
!X perform different tasks. If an error occurs on one process it will  
!X be broadcast to all others before the error is propagated  
!X upwards. Replace RAISE_ERROR with BCAST_RAISE_ERROR and PASS_ERROR  
!X with BCAST_PASS_ERROR. Additionally, BCAST_CHECK_ERROR must be  
!X called on the processes in which no error has occured.  See  
!X CInOutput read() for an example usage of these macros.  
!X  
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
  
#define BCAST_ASSERT(condition, message, error, mpi) if (.not. (condition)) then ; BCAST_RAISE_ERROR(message, error, mpi) ; endif  
  
#define BCAST_RAISE_ERROR(message, error, mpi) if (present(error)) call bcast(mpi, error); RAISE_ERROR(message, error)  
  
#define BCAST_RAISE_ERROR_WITH_KIND(kind, message, error, mpi) if (present(error)) then; error = kind; call bcast(mpi, error); end if; RAISE_ERROR_WITH_KIND(kind, message, error)  
  
#define BCAST_PASS_ERROR(error, mpi) if (present(error)) then; if (error /= ERROR_NONE) call bcast(mpi, error); endif; PASS_ERROR(error)  
  
#define BCAST_CHECK_ERROR(error, mpi) if (present(error)) then; call bcast(mpi, error); if (error /= ERROR_NONE) then; RAISE_ERROR_WITH_KIND(error, "An error occured on another MPI process", error); endif; endif  
  
  
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
!X  
!X Delayed errors - for OpenMP loops  
!X  
!X A subroutine currently in an OpenMP section cannot be quit using  
!X the *return* statement. Hence, the error flag is set using  
!X RAISE_DELAYED_ERROR and TRACE_DELAYED_ERROR. After the OpenMP section  
!X has finished, INVOKE_DELAYED_ERROR will raise the error and exit  
!X the current subroutine if an error occured in the OpenMP section.  
!X  
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
  
#define RAISE_DELAYED_ERROR(message, error_loc)   if (error_loc == ERROR_NONE) then ; call push_error_with_info(message, __FILE__, __LINE__) ; error_loc = ERROR_UNSPECIFIED ; endif  
  
#define TRACE_DELAYED_ERROR(error_loc)  if (error_loc /= ERROR_NONE) then ; call push_error(__FILE__, __LINE__) ; endif  
  
#define TRACE_DELAYED_ERROR_WITH_INFO(message, error_loc)  if (error_loc /= ERROR_NONE) then ; call push_error_with_info(message, __FILE__, __LINE__) ; endif  
  
#define INVOKE_DELAYED_ERROR(error_loc, error)  if (error_loc /= ERROR_NONE) then ; call push_error(__FILE__, __LINE__) ; if (present(error)) then ; error = error_loc ; else ; call error_abort(error) ; endif ; endif  
  
