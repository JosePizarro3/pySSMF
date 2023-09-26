!*****************************************************************************************
!*** Copyright: Dr. José M. Pizarro.
!***
!*** Licensed under the Apache License, Version 2.0 (the "License");
!*** you may not use this file except in compliance with the License.
!*** You may obtain a copy of the License at
!***
!***     http://www.apache.org/licenses/LICENSE-2.0
!***
!*** Unless required by applicable law or agreed to in writing, software
!*** distributed under the License is distributed on an "AS IS" BASIS,
!*** WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!*** See the License for the specific language governing permissions and
!*** limitations under the License.
!*****************************************************************************************


! Subroutine to diagonalize a complex Hermitian matrix.
! Inputs:
!   n_k: Number of k points
!   n_orb: Number of orbitals
!   hamiltonian: Input Hamiltonian (complex Hermitian matrix of shape (n_k, n_orb, n_orb))
!
! Outputs:
!   eigenvectors_out: Eigenvectors (complex Hermitian matrix of shape (n_k, n_orb, n_orb))
!   eigenvalues_out: Eigenvalues (real matrix of shape (n_k, n_orb))
!
! Note: This subroutine uses LAPACK's zheev routine for diagonalization.

program diagonalize
    implicit none

    ! Input dimensions
    integer, parameter :: max_k = 300   ! number of k points
    integer, parameter :: max_orb = 10 ! number of orbitals
    integer :: n_k, n_orb
    ! Input Hamiltonian (complex Hermitian matrix)
    complex*8 :: hamiltonian(max_k, max_orb, max_orb)

    ! Output eigensystem
    complex*16, allocatable :: eigenvectors_out(:, :, :)
    real*8, allocatable :: eigenvalues_out(:, :)

    ! General integers for loops
    integer k, i_orb, j_orb, i_band

    ! LAPACK parameters
    complex*16, allocatable :: eigenvectors(:, :)
    real*8, allocatable :: eigenvalues(:)
    integer :: lwork
    real*8, allocatable :: rwork(:)
    complex*16, allocatable :: work(:)
    character jobz, uplo
    integer lda, info
    jobz = 'V'
    uplo = 'U'

    ! Read input from stdin
    read(*, *) n_k, n_orb
    do k = 1, n_k
        do i_orb = 1, n_orb
            do j_orb = 1, n_orb
                read(*, *) hamiltonian(k, i_orb, j_orb)
            end do
        end do
    end do

    ! Allocate memory for output
    allocate(eigenvectors_out(n_k, n_orb, n_orb))
    allocate(eigenvalues_out(n_k, n_orb))
    allocate(eigenvectors(n_orb, n_orb))
    allocate(eigenvalues(n_orb))
    allocate(rwork(3 * n_orb - 2))
    allocate(work(2 * n_orb - 1))
    lda = n_orb
    lwork = 2 * n_orb - 1

    ! Diagonalize for each n_k
    do k = 1, n_k
        do i_orb = 1, n_orb
            do j_orb = 1, n_orb
                eigenvectors(i_orb, j_orb) = hamiltonian(k, i_orb, j_orb)
            end do
        end do
        call zheev(jobz, uplo, lda, eigenvectors, lda, eigenvalues, work, lwork, rwork, info)

        ! Check for errors
        if (info /= 0) then
            write(*, *) 'Error: zheev diagonalization routine failed with info = ', info
            stop
        end if

        ! Store results in eigenvectors_out and eigenvalues_out
        do i_band = 1, n_orb
            do i_orb = 1, n_orb
                eigenvectors_out(k, i_orb, i_band) = eigenvectors(i_orb, i_band)
            end do
        end do
    end do

    ! Write output to stdout
    write(*, *) n_k, n_orb
    do k = 1, n_k
        do i_band = 1, n_orb
            write(*, *) eigenvalues_out(k, i_band)
        end do
    end do

    ! Deallocate memory
    deallocate(eigenvectors_out)
    deallocate(eigenvalues_out)
end program diagonalize
