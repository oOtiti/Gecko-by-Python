# ================================================================
# the subroutine sort_string(string) return the array "string" provided
# as input as a sorted array (alphabetical order).
# ================================================================

# 排序阈值：快速排序处理到该大小后，切换为插入排序（原Fortran参数QSORT_THRESHOLD）
QSORT_THRESHOLD = 16
# 栈大小：处理数组索引的最大深度（32位索引足够覆盖常规场景，原Fortran参数QSORT_STACK_SIZE）
QSORT_STACK_SIZE = 32

def less_than(a: int, b: int, string: list[str]) -> bool:
    """
    LOGICAL FUNCTION less_than(a,b,string)
      INTEGER, INTENT(IN) :: a, b
      CHARACTER(LEN=*), INTENT(IN) :: string(:)
      less_than = (string(a) < string(b))
    END FUNCTION
    """
    return string[a - 1] < string[b - 1]

def swap(a: int, b: int, string: list[str]) -> None:
    """
    SUBROUTINE swap(a,b,string)
      INTEGER, INTENT(IN) :: a, b
      CHARACTER(LEN=*),INTENT(INOUT) :: string(:)
      CHARACTER(LEN=LEN(string(1)))  :: hold_string
      hold_string = string(a)
      string(a) = string(b)
      string(b) = hold_string
      RETURN
    END SUBROUTINE
    """
    hold_string = string[a - 1]
    string[a - 1] = string[b - 1]
    string[b - 1] = hold_string

def r_shift(left: int, right: int, string: list[str]) -> None:
    """
    SUBROUTINE r_shift(left,right,string)
      INTEGER, INTENT(IN) :: left, right
      CHARACTER(LEN=*),INTENT(INOUT) :: string(:)
      INTEGER  ::  i
      CHARACTER(LEN=LEN(string(1)))  :: hold_string
      hold_string = string(right)
      DO i=right, left+1, -1
        string(i) = string(i-1)
      ENDDO
      string(left) = hold_string
    END SUBROUTINE
    """
    if left >= right:
        return
    hold_string = string[right - 1]
    # 从right到left+1反向循环（1-based转0-based）
    for i in range(right, left, -1):
        string[i - 1] = string[i - 2]
    string[left - 1] = hold_string

def sort_string(string: list[str]) -> None:
    """
    SUBROUTINE sort_string(string)
!======================================================================
! Fast in-line QSORT+INSERTION SORT for Fortran.
! Author: Joseph M. Krahn
! FILE: qsort_inline.inc
! PURPOSE:
! Generate a custom array sort procedure for a specific type,
! without the comparison-callback overhead of a generic sort procedure.
! This is essentially the same as an in-line optimization, which generally
! is not feasible for a library-based generic sort procedure.
!
! This implementation is as generic as possible, while avoiding the need
! for a code pre-processor. The success of this approach assumes that
! internal procedures are always in-lined with optimized Fortran compilation.
!
! USAGE:
! This file contains the sort subroutine body. You must supply
! an integer parameter QSORT_THRESHOLD, and internal procedures:
!    subroutine INIT()
!    logical function LESS_THAN(a,b)
!    subroutine SWAP(a,b)
!    subroutine RSHIFT(left,right)
!
! Descriptions:
!
! SUBROUTINE INIT()
!   Any user initialization code. This is needed because executable
!   statements cannot precede this code, which begins with declarations.
!   In many cases, this is just an empty procedure.
!
! LOGICAL FUNCTION LESS_THAN(a,b)
!   Return TRUE if array member 'a' is less than array member 'b'.
!   Only a TRUE value causes a change in sort order. This minimizes data
!   manipulation, and maintains the original order for values that are
!   equivalent by the sort comparison. It also avoids the need to
!   distinguish equality from greater-than.
!
! SUBROUTINE SWAP(A,B)
!   Swap array members 'a' and 'b'
!
! SUBROUTINE RSHIFT(LEFT,RIGHT)
!   Perform a circular shift of array members LEFT through RIGHT,
!   shifting the element at RIGHT back to the position at LEFT.
!
! QSORT_THRESHOLD:
!   The QSORT is used down to the QSORT_THRESHOLD size sorted blocks.
!   Then insertion sort is used for the remainder, because it is faster
!   for small sort ranges. The optimal size is not critical. Most of
!   the benefit is in blocks of 8 or less, and values of 16 to 128
!   are generally about equal speed. However, the optimal value
!   depends a lot on the hardware and the data being sorted, so this
!   is left as a tunable parameter for cases where ther is an
!   effect on performance.
!
!---------------------------------------------------------------------
! NOTES:
! The procedure uses a optimized combination of QSORT and INSERTION
! sorting. The algorithm is based on code used in GLIBC.
! A stack is used in place of recursive calls. The stack size must
! be at least as big as the number of bits in the largest array index.
!
! Sorting vectors of a multidimensional allocatable array can be
! VERY slow. In this case, or with large derived types, it is better
! to sort a simple derived type of key/index pairs, then reorder
! tha actual data using the sorted indices.
!---------------------------------------------------------------------
    """
    array_size = len(string)
    if array_size <= 1:
        return

    # A stack of 32 can handle the entire extent of a 32-bit
    # index, so this value is fixed. If you have 64-bit indexed
    # arrays, which might contain more thant 2^32 elements, this
    # should be set to 64.
    stack = [(0, 0)] * QSORT_STACK_SIZE  # 每个元素存储(low, high)，模拟原Fortran的qsort_stack类型
    stack_top = 0
    right_size = 0
    left_size = 0
    mid = 0
    left = 0
    right = 0
    low = 0
    high = 0

    if array_size > QSORT_THRESHOLD:
        low = 1
        high = array_size
        stack_top = 0

        # QSORT_LOOP:
        while True:
            mid = (low + high) // 2
            # if (LESS_THAN (mid, low)) then
            if less_than(mid, low, string):
                # call SWAP(mid,low)
                swap(mid, low, string)
            # if (LESS_THAN (high, mid)) then
            if less_than(high, mid, string):
                # call SWAP(high,mid)
                swap(high, mid, string)
                # if (LESS_THAN (mid, low)) then
                if less_than(mid, low, string):
                    # call SWAP(mid,low)
                    swap(mid, low, string)
            left = low + 1
            right = high - 1

            # COLLAPSE_WALLS:
            while True:
                # do while (LESS_THAN (left, mid))
                while less_than(left, mid, string):
                    left += 1
                # do while (LESS_THAN (mid, right))
                while less_than(mid, right, string):
                    right -= 1
                if left < right:
                    # call SWAP(left,right)
                    swap(left, right, string)
                    if mid == left:
                        mid = right
                    elif mid == right:
                        mid = left
                    left += 1
                    right -= 1
                else:
                    if left == right:
                        left += 1
                        right -= 1
                    break  # exit COLLAPSE_WALLS

            # Set up indices for the next iteration.
            # Determine left and right partition sizes.
            # Defer partitions smaller than the QSORT_THRESHOLD.
            # If both partitions are significant,
            # push the larger one onto the stack.
            right_size = right - low
            left_size = high - left
            if right_size <= QSORT_THRESHOLD:
                if left_size <= QSORT_THRESHOLD:
                    # Ignore both small partitions: Pop a partition or exit.
                    if stack_top < 1:
                        break  # exit QSORT_LOOP
                    low, high = stack[stack_top - 1]
                    stack_top -= 1
                else:
                    # Ignore small left partition.
                    low = left
            elif left_size <= QSORT_THRESHOLD:
                # Ignore small right partition.
                high = right
            elif right_size > left_size:
                # Push larger left partition indices.
                stack_top += 1
                if stack_top > QSORT_STACK_SIZE:
                    raise OverflowError("QSORT stack overflow: increase QSORT_STACK_SIZE")
                stack[stack_top - 1] = (low, right)
                low = left
            else:
                # Push larger right partition indices.
                stack_top += 1
                if stack_top > QSORT_STACK_SIZE:
                    raise OverflowError("QSORT stack overflow: increase QSORT_STACK_SIZE")
                stack[stack_top - 1] = (left, high)
                high = right

    # Sort the remaining small partitions using insertion sort, which should
    # be faster for partitions smaller than the appropriate QSORT_THRESHOLD.

    # First, find smallest element in first QSORT_THRESHOLD and place it at
    # the array's beginning. This places a lower bound 'guard' position, and
    # speeds up the inner loop below, because it will not need a lower-bound
    # test.
    low = 1
    high = array_size

    # left is the MIN_LOC index here:
    left = low
    max_right = min(low + QSORT_THRESHOLD, high)
    for right_val in range(low + 1, max_right + 1):
        # if (LESS_THAN(right_val, left)) left=right_val
        if less_than(right_val, left, string):
            left = right_val
    # if (left/=low) call SWAP(left,low)
    if left != low:
        swap(left, low, string)

    # Insertion sort, from left to right.
    # (assuming that the left is the lowest numbered index)
    # insertion_sort:
    for right_val in range(low + 2, high + 1):
        left = right_val - 1
        # if (LESS_THAN(right_val, left)) then
        if less_than(right_val, left, string):
            # do while (LESS_THAN(right_val, left-1))
            while less_than(right_val, left - 1, string):
                left -= 1
            # call R_SHIFT(left, right_val)
            r_shift(left, right_val, string)