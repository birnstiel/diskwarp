! inspired by this: https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
! and this https://courses.cs.washington.edu/courses/csep557/10au/lectures/triangle_intersection.pdf
! locating the points in a skewed/projected grid that can have overlap.


MODULE fmodule
   IMPLICIT NONE


CONTAINS
! -------------------------------------------------------------------
   FUNCTION nan()
      DOUBLE PRECISION nan
      nan = 0.0
      nan = 0.0 / nan
   END FUNCTION

! -------------------------------------------------------------------
   FUNCTION signof(p1, p2, p3)
      IMPLICIT NONE
      DOUBLE PRECISION :: signof
      DOUBLE PRECISION, dimension(2) :: p1, p2, p3
      signof = (p1(1) - p3(1))*(p2(2) - p3(2)) - (p2(1) - p3(1))*(p1(2) - p3(2))
      return
   END FUNCTION

! -------------------------------------------------------------------
!checks if 2D-point pt is inside triangle spanned by points v1, v2, v3
   FUNCTION PointInTriangle(pt, v1, v2, v3)
      logical :: PointInTriangle
      DOUBLE PRECISION :: d1, d2, d3
      DOUBLE PRECISION, dimension(2) :: v1, v2, v3, pt
      logical :: has_neg, has_pos

      if (ANY(pt .ne. pt)) then
         PointInTriangle = .false.
         RETURN
      end if

      d1 = signof(pt, v1, v2)
      d2 = signof(pt, v2, v3)
      d3 = signof(pt, v3, v1)

      has_neg = (d1 .lt. 0) .or. (d2 .lt. 0) .or. (d3 .lt. 0)
      has_pos = (d1 .gt. 0) .or. (d2 .gt. 0) .or. (d3 .gt. 0)

      PointInTriangle = .not. (has_neg .and. has_pos)
   END FUNCTION

! -------------------------------------------------------------------
! same for square: will separate into two triangles and check if point pt is in one of them
   logical FUNCTION PointInSquare(pt, v1, v2, v3, v4)
   DOUBLE PRECISION, dimension(2) :: v1, v2, v3, v4, pt

      PointInSquare = PointInTriangle(pt, v1, v2, v3) .or. PointInTriangle(pt, v1, v3, v4)
   END FUNCTION
! -------------------------------------------------------------------
! find the grid cells that are intersected and pick the front one
   SUBROUTINE findpoint(xi, yi, pt, nx, ny, matches, i_match)
      IMPLICIT NONE
      INTEGER, INTENT(in):: nx, ny
      DOUBLE PRECISION, INTENT(in), dimension(nx, ny) :: xi, yi
      DOUBLE PRECISION, INTENT(in), dimension(2) :: pt

      INTEGER, INTENT(out) :: i_match
      INTEGER, PARAMETER :: n_match = 10 ! there should just be a hand full of matches
      INTEGER, INTENT(out), DIMENSION(n_match, 2) :: matches

      DOUBLE PRECISION, dimension(2) :: v1, v2, v3, v4
      INTEGER :: ix, iy

      i_match = 0
      ! loop over all the desired pixels and find which 4-sided polygon includes the point
      do ix = 1, nx - 1
         do iy = 1, ny - 1
            v1(1) = xi(ix, iy)
            v1(2) = yi(ix, iy)
            v2(1) = xi(ix + 1, iy)
            v2(2) = yi(ix + 1, iy)
            v3(1) = xi(ix + 1, iy + 1)
            v3(2) = yi(ix + 1, iy + 1)
            v4(1) = xi(ix, iy + 1)
            v4(2) = yi(ix, iy + 1)

            if (PointInSquare(pt, v1, v2, v3, v4)) then
               i_match = i_match + 1
               matches(i_match, 1) = ix
               matches(i_match, 2) = iy
            END if
         end do
      end do

   END SUBROUTINE

! -------------------------------------------------------------------
   FUNCTION cross(a, b)
      DOUBLE PRECISION, DIMENSION(3) :: cross
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: a, b
    
      cross(1) = a(2) * b(3) - a(3) * b(2)
      cross(2) = a(3) * b(1) - a(1) * b(3)
      cross(3) = a(1) * b(2) - a(2) * b(1)
    END FUNCTION cross

! -------------------------------------------------------------------
    ! this finds the intersection of a z-directed ray through the
    ! point P with a plane spanned by points A, B, C.
    FUNCTION intersect_triangle(A, B, C, P)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3) :: intersect_triangle
      DOUBLE PRECISION, INTENT(IN) :: A(3), B(3), C(3), P(3)
      DOUBLE PRECISION :: n(3), ray(3), n_dot_ray, d_plane, t

      ! normal of the supporting plane
      n = cross((B - A), (C - A))
      n = n / norm2(n)
      d_plane = DOT_PRODUCT(n, A)

      ray = (/ 0, 0, 1 /)
      n_dot_ray = DOT_PRODUCT(n, ray)
      if (n_dot_ray .eq. 0.0) then
         ! parallel to plane: no intersection
         intersect_triangle = intersect_triangle * nan()
      else
         t = (d_plane - DOT_PRODUCT(n, P)) / n_dot_ray
         intersect_triangle = P + t * ray
      endif

   END FUNCTION

! -------------------------------------------------------------------
! subroutine to bary-center interpolate at the point Q on a triangle
! spanned by points (A, B, C) with their values (VA, VB, VC).
! nq quantities can be interpolated, so VA, VB, VC have shape (nq)
   subroutine barycenter_interpolate(A, B, C, Q, VA, VB, VC, result, nq)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nq
      DOUBLE PRECISION, INTENT(IN) :: A(3), B(3), C(3), Q(3), VA(nq), VB(nq), VC(nq)
      DOUBLE PRECISION, INTENT(OUT) :: result(nq)
      DOUBLE PRECISION :: alpha, beta, gamma, n(3)
      ! normal of the supporting plane
      n  = cross((B-A), (C-A))
      n  = n / norm2(n)

      alpha = DOT_PRODUCT(cross(C-B, Q-B), n) / DOT_PRODUCT(cross(B-A, C-A), n)
      beta =  DOT_PRODUCT(cross(A-C, Q-C), n) / DOT_PRODUCT(cross(B-A, C-A), n)
      gamma =  DOT_PRODUCT(cross(B-A, Q-A), n) / DOT_PRODUCT(cross(B-A, C-A), n)

      result = alpha * VA + beta * VB + gamma * VC
   END subroutine

! -------------------------------------------------------------------
! find the grid cells that are intersected and pick the front one
   SUBROUTINE intersect_surface(xi, yi, zi, vi, pt, nx, ny, matches, n_match, n_max, nq)
      IMPLICIT NONE
      INTEGER, INTENT(IN):: nx, ny, n_max, nq
      DOUBLE PRECISION, INTENT(in), dimension(nx, ny) :: xi, yi, zi
      DOUBLE PRECISION, INTENT(in), dimension(nx, ny, nq) :: vi
      DOUBLE PRECISION, INTENT(in), dimension(2) :: pt

      DOUBLE PRECISION, INTENT(OUT), DIMENSION(n_max, 3+nq) :: matches
      INTEGER,  INTENT(OUT) :: n_match
      DOUBLE PRECISION :: Q(3), V(nq)

      DOUBLE PRECISION, dimension(3) :: p1, p2, p3, p4
      DOUBLE PRECISION :: v1(nq), v2(nq), v3(nq), v4(nq), p(3)
      INTEGER :: ix, iy
      
      LOGICAL :: match

      p(1:2) = pt
      p(3) = 0.0
      n_match = 0

      ! loop over all the square of the grid and look for intersections
      do iy = 1, ny - 1
         do ix = 1, nx - 1
            p1(1) = xi(ix, iy)
            p1(2) = yi(ix, iy)
            
            p2(1) = xi(ix + 1, iy)
            p2(2) = yi(ix + 1, iy)
            
            p3(1) = xi(ix + 1, iy + 1)
            p3(2) = yi(ix + 1, iy + 1)

            p4(1) = xi(ix, iy + 1)
            p4(2) = yi(ix, iy + 1)

            match = .false.

            ! if we are clearly left or right of the cell, skip to next
            if      (pt(1) .gt. MAX(p1(1), p2(1), p3(1), p4(1))) then
               cycle
            else if (pt(1) .lt. MIN(p1(1), p2(1), p3(1), p4(1))) then
               cycle
            else if (pt(2) .lt. MIN(p1(2), p2(2), p3(2), p4(2))) then
               cycle
            else if (pt(2) .gt. MAX(p1(2), p2(2), p3(2), p4(2))) then
               cycle
            endif

            ! with that test out of the way, we can get the other relevant quantities ...
            p1(3) = zi(ix, iy)
            p2(3) = zi(ix + 1, iy)
            p3(3) = zi(ix + 1, iy + 1)
            p4(3) = zi(ix, iy + 1)
            
            v1 = vi(ix, iy, :)
            v2 = vi(ix + 1, iy, :)
            v3 = vi(ix + 1, iy + 1, :)
            v4 = vi(ix, iy + 1, :)

            ! and check for intersections in detail

            if (PointInTriangle((/pt(1),pt(2)/), (/p1(1),p1(2)/), (/p2(1),p2(2)/), (/p3(1),p3(2)/))) then
               match = .true.
               ! interpolate the point on the triangle and return its value and z-coordinate
               Q = intersect_triangle(p1, p2, p3, p)
               call barycenter_interpolate(p1, p2, p3, Q, v1, v2, v3, V, nq)
               
            else if (PointInTriangle((/pt(1),pt(2)/), (/p1(1),p1(2)/), (/p3(1),p3(2)/), (/p4(1),p4(2)/))) then
               match = .true.
               ! interpolate the point on the triangle and return its value and z-coordinate
               Q = intersect_triangle(p1, p3, p4, p)
               call barycenter_interpolate(p1, p3, p4, Q, v1, v3, v4, V, nq)
            ENDIF
            if (match) then
               n_match = n_match + 1
               matches(n_match, 1:3) = Q
               matches(n_match, 4:3+nq) = V
            ENDIF
         end do
      end do

   END SUBROUTINE

! -------------------------------------------------------------------
! function to project the 3d data along zi:
! xi, yi, zi : 2D matrices of shape (nx, ny) giving the x, y, z positions and values at those positions
! vi: 3D matrix of the quantities that should be interpolated, shape (nx, ny, nq) for nq quantities.
! img_x, img_y : 2D matrices of shape (nix, niy) giving the x, y, z positions of the pixels where the projections should be carried out
! returns:
! - img_z : 2D matrix of shape (nix, niy) giving the intersection z-value at img_x, img_y
! - img_v : 2D matrix of shape (nix, niy, nq) giving the intersection value at img_x, img_y, for all nq quantities
   SUBROUTINE interpolate_grid(xi, yi, zi, vi, img_x, img_y, img_z, img_v, nix, niy, nx, ny, nq)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nix, niy, nx, ny, nq
      DOUBLE PRECISION, INTENT(IN) :: img_x(nix, niy), img_y(nix, niy), xi(nx, ny), yi(nx, ny), zi(nx, ny), vi(nx, ny, nq)
      DOUBLE PRECISION, INTENT(OUT) :: img_z(nix, niy), img_v(nix, niy, nq)
      INTEGER, PARAMETER :: n_max = 10
      DOUBLE PRECISION :: matches(n_max, 3 + nq)
      INTEGER :: numbers(n_max)
      DOUBLE PRECISION :: x, y, z_vals(n_max)
      INTEGER :: ix, iy, loc(1), n_match, i

      img_z = 0.0
      img_v = 0.0

      do i = 1, n_max
         numbers(i) = i
      end do

      do ix = 1, nix
         do iy = 1, niy
            x = img_x(ix, iy)
            y = img_y(ix, iy)

            ! for each image pixel, we find the intersection
            call intersect_surface(xi, yi, zi, vi, (/x, y/), nx, ny, matches, n_match, n_max, nq)

            if (n_match .eq. 0) then
               img_z(ix, iy) = nan()
               img_v(ix, iy, :) = nan()
               continue
            endif

            z_vals = matches(:, 3)

            ! we find index where the z-values (height along point of view) is largest
            loc = MAXLOC(z_vals, DIM=1, MASK=n_match .ge. numbers)

            ! get the coordinates
            img_z(ix, iy) = matches(loc(1), 3)
            img_v(ix, iy, :) = matches(loc(1), 4:nq+3)

         end do
      end do

   end SUBROUTINE

! -------------------------------------------------------------------

   subroutine apply_matrix(p, warp, twist, pout)

      DOUBLE PRECISION, INTENT(in) :: p(3), warp
      DOUBLE PRECISION, INTENT(IN), optional :: twist
      DOUBLE PRECISION, INTENT(OUT) :: pout(3)
      ! the default values
      DOUBLE PRECISION :: phi_t = 0d0
      DOUBLE PRECISION :: x, y, z
      DOUBLE PRECISION :: cosw, sinw, cost, sint
      if(present(twist)) phi_t=twist
   
      x = p(1)
      y = p(2)
      z = p(3)
   
      cosw = cos(warp)
      sinw = sin(warp)
      
      cost = cos(phi_t)
      sint = sin(phi_t)
      
      !Rx(warp) Rz(-phi_t) Goes in the correct z direction
      pout(1) = x*cost + y*sint
      
      pout(2) = -x*sint*cosw + y*cost*cosw - z*sinw
   
      pout(3) = -x*sint*sinw + y*sinw*cost + z*cosw
   
   end subroutine
   
   ! -------------------------------------------------------------------
   
   subroutine apply_matrix2D(p, warp, twist, pout, nr, nphi)
   
   INTEGER, INTENT(in) :: nr, nphi
   DOUBLE PRECISION, INTENT(in) :: p(nr, nphi, 3), warp(nr)
   DOUBLE PRECISION, INTENT(IN), optional :: twist(nr)
   DOUBLE PRECISION, INTENT(OUT) :: pout(nr, nphi, 3)
   ! the default values
   DOUBLE PRECISION :: phi_t(nr)
   integer :: ir, iphi
   
   if(present(twist)) then 
      phi_t=twist
   else
      phi_t=0d0
   end if
   
   do ir = 1, nr
      do iphi = 1, nphi
         call apply_matrix(p(ir, iphi, :), warp(ir), phi_t(ir), pout(ir, iphi, :))
      end do
   end do
   
   
   end subroutine
   
   ! -------------------------------------------------------------------
   

subroutine test_module(xi, yi, zi, vi, img_x, img_y, img_z, img_v, nq)
   IMPLICIT NONE
   INTEGER, parameter :: nix=30, niy=20
   INTEGER, INTENT(IN) :: nq
   DOUBLE PRECISION, INTENT(OUT) :: xi(2, 2), yi(2, 2), zi(2, 2), vi(2, 2, nq)
   DOUBLE PRECISION, INTENT(OUT) :: img_x(nix, niy), img_y(nix, niy), img_z(nix, niy), img_v(nix, niy, nq)
   INTEGER:: ix, iy

   xi(1, 1) = 1.0
   xi(2, 1) = 1.5
   xi(2, 2) = 3.5
   xi(1, 2) = 2.0

   yi(1, 1) = 1.0
   yi(2, 1) = 3.0
   yi(2, 2) = 2.0
   yi(1, 2) = 0.0

   zi(1, 1) = 1.0
   zi(2, 1) = 1.0
   zi(2, 2) = 2.0
   zi(1, 2) = 2.0

   vi(1, 1, :) = 1.0
   vi(2, 1, :) = 2.0
   vi(2, 2, :) = 3.0
   vi(1, 2, :) = 2.0

   do ix = 1, nix
      do iy = 1, niy
         img_x(ix, iy) = (ix - 1.0) / (nix - 1.0) * 4.0
         img_y(ix, iy) = (iy - 1.0) / (niy - 1.0) * 4.0
      enddo
   enddo

   call interpolate_grid(xi, yi, zi, vi, img_x, img_y, img_z, img_v, nix, niy, 2, 2, nq)

end subroutine 


END MODULE fmodule