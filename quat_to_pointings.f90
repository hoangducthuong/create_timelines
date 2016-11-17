! Compile this with
!
!     f2py3 -c -m quat_to_pointings --f90flags=-std=f2003 quat_to_pointings.f90
!

subroutine quaternion_to_axis_angle(q, axis, angle)

  implicit none

  real(kind=8), dimension(:, :), intent(in)           :: q
  real(kind=8), dimension(size(q, 1), 3), intent(out) :: axis
  real(kind=8), dimension(size(q, 1)), intent(out)    :: angle

  integer :: idx

  do idx = 1, size(q, 1)
     angle(idx) = 2 * acos(q(idx, 4))
     if (angle(idx) == 0.0) then
        axis(idx, :) = 0.0
     else
        axis(idx, :) = q(idx, 1:3) / sin(angle(idx)/2.0)
     endif
  end do

end subroutine quaternion_to_axis_angle
