module typedefs

  type dmat
    integer :: n
    real(8), pointer :: h(:, :)
  end type dmat

end module typedefs
