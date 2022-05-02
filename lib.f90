function trial(n, trials) result (successes)
  implicit none

  ! constants
  real, parameter, dimension(25) :: probs = (/ &
    0.1202, 0.0910, 0.0812, 0.0149, 0.0271, &
    0.0432, 0.0230, 0.0203, 0.0592, 0.0731, &
    0.0010, 0.0069, 0.0398, 0.0261, 0.0695, &
    0.0768, 0.0182, 0.0011, 0.0602, 0.0628, &
    0.0288, 0.0111, 0.0209, 0.0017, 0.0211 /)

  integer, intent(in) :: n
  integer, intent(in) :: trials
  integer :: successes
  real, dimension(n) :: rands
  integer, dimension(26) :: counts

  ! Looping variables
  integer :: t
  integer :: l
  integer :: i

  successes = 0

  trialloop: do t = 1, trials
    call random_number(rands)
    counts = (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
    do l = 1, n
      do i = 1, 25
        rands(l) = rands(l) - probs(i)
        if (rands(l) <= 0) then
          counts(i) = counts(i) + 1
          exit
        else if (i == 25) then
          counts(26) = counts(26) + 1
        end if
      end do
    end do

    if (counts(2) >= counts(1)) then
      cycle trialloop
    end if
    do i = 3, 26
      if (counts(i) >= counts(2)) then
        cycle trialloop
      end if
    end do

    successes = successes + 1
  end do trialloop
end function trial