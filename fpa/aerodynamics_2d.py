"""
This module enables XLFR5 data to numpy.array.
"""

import numpy as np
import os
# import matplotlib.pyplot as plt


def read_xflr5_data(wing_foil_name):
    reynolds_numbers = np.arange(0, 1.01, 0.1)
    headers = ["alpha", "CL", "CD", "CDp", "Cm", "Top Xtr", "Bot Xtr", "Cpmin", "Chinge", "XCp"]
    tmpdict = {i: [] for i in headers}

    for reynolds_number in reynolds_numbers:
        data = np.genfromtxt(os.path.join(wing_foil_name, "{:.2f}.txt".format(reynolds_number)),
                             skip_header=11).transpose()
        for header, datum in zip(headers, data):
            tmpdict[header].append(datum)
    return {i: np.array(tmpdict[i]) for i in tmpdict}
# print tmpdict["Cm"]

# for j in range(len(reynolds_numbers)):
#     plt.plot(tmpdict["alpha"][j], tmpdict["Cm"][j], label="Reynolds={}".format(j))
# plt.legend(loc=2, fontsize="x-small")
# plt.show()


if __name__ == '__main__':
    pass

