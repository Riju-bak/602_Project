def find_ind_val(arr, val):
    for i in range(len(arr)):
        if abs(arr[i] - val) < 0.01:
            return i
    return -1
