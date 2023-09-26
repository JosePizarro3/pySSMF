

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
        do i_orb = 1, n_orb
            do i_band = 1, n_orb
                write(*, *) eigenvectors_out(k, i_orb, i_band)
            end do
        end do
        do i_band = 1, n_orb
            write(*, *) eigenvalues_out(k, i_band)
        end do
    end do