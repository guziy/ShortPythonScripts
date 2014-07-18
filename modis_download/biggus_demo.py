__author__ = 'huziy'

import numpy as np
import biggus

class ArrForDate:
    def __init__(self):
        self.shape = (100, 100)
        self.dtype = np.dtype(np.float32)
        self.fill_value = -1

    def __getitem__(self, item):
        return np.ones(self.shape, dtype=self.dtype)[item]


def main():
    arr_stack = biggus.ArrayStack(
        np.array([biggus.NumpyArrayAdapter(ArrForDate()) for _ in range(20)])
    )

    print arr_stack.shape
    import matplotlib.pyplot as plt
    plt.pcolormesh(biggus.mean(arr_stack, axis=0).ndarray())
    plt.show()



if __name__ == '__main__':
    main()
