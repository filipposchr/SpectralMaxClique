import main as m

#n=1000, a=1/3, k=1,..,6
    #k=1
probabilities = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.13, 0.15, 0.18, 0.2]
m.run_tests(1000, 1/3, 1, probabilities)
    #k=2
probabilities = [0.06, 0.07, 0.08, 0.1, 0.11, 0.12, 0.14, 0.16, 0.18, 0.2]
m.run_tests(1000, 1/3, 2, probabilities)
    #k=3
probabilities = [0.06, 0.07, 0.08, 0.1, 0.11, 0.12,  0.15, 0.2, 0.21, 0.22, 0.24, 0.26]
m.run_tests(1000, 1 / 3, 3, probabilities)
    #k=4
probabilities = [0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.21, 0.22, 0.23, 0.25, 0.28]
m.run_tests(1000, 1 / 3, 4, probabilities)
    #k=5
probabilities = [0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.21, 0.22, 0.23, 0.25, 0.26, 0.28, 0.3, 0.32]
m.run_tests(1000, 1 / 3, 5, probabilities)
    #k=6
probabilities = [0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28,  0.29, 0.3, 0.32, 0.34, 0.36]
m.run_tests(1000, 1 / 3, 6, probabilities)


#n=1000, a=2/3, k=1,..,7
    #k=1
probabilities = [0.003, 0.005, 0.01, 0.02, 0.03, 0.035, 0.04, 0.05, 0.06]
m.run_tests(1000, 2/3, 1, probabilities)
    #k=2
probabilities = [0.005, 0.006, 0.007, 0.008, 0.01, 0.015, 0.02, 0.03, 0.035, 0.04, 0.045, 0.05, 0.06, 0.07]
m.run_tests(1000, 2/3, 2, probabilities)
    #k=3
probabilities = [0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08]
m.run_tests(1000, 2/3, 3, probabilities)
    #k=4
probabilities = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.09, 0.1]
m.run_tests(1000, 2/3, 4, probabilities)
    #k=5
probabilities = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.075, 0.08, 0.09, 0.1]
m.run_tests(1000, 2/3, 5, probabilities)
    #k=6
probabilities = [0.03, 0.035, 0.04,  0.042, 0.044, 0.046, 0.05, 0.06, 0.07, 0.074, 0.076, 0.078, 0.08, 0.082, 0.084]
m.run_tests(1000, 2/3, 6, probabilities)
    #k=7
probabilities = [0.04, 0.042, 0.043, 0.044, 0.046, 0.048, 0.05, 0.06, 0.07, 0.075, 0.08, 0.081, 0.082, 0.084, 0.086, 0.088]
m.run_tests(1000, 2/3, 7, probabilities)


#n=3000, a=1, k=1,..,8
    #k=1
probabilities = [0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.0012, 0.0013, 0.0014, 0.0015, 0.0016, 0.0018, 0.002]
m.run_tests(3000, 1, 1, probabilities)
    #k=2
probabilities = [0.001, 0.0015, 0.002, 0.0025, 0.003, 0.004, 0.0045, 0.005, 0.006]
m.run_tests(3000, 1, 2, probabilities)
    #k=3
probabilities = [0.003, 0.0035, 0.004, 0.005, 0.006, 0.0065, 0.007, 0.0075, 0.008, 0.009]
m.run_tests(3000, 1, 3, probabilities)
    #k=4
probabilities = [0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.012, 0.014]
m.run_tests(3000, 1, 4, probabilities)
    #k=5
probabilities = [0.005, 0.006, 0.007, 0.008,  0.009, 0.01, 0.011, 0.012, 0.014, 0.016]
m.run_tests(3000, 1, 5, probabilities)
    #k=6
probabilities = [0.007, 0.008, 0.009, 0.01, 0.011, 0.012,  0.013, 0.014, 0.016, 0.018]
m.run_tests(3000, 1, 6, probabilities)
    #k=7
probabilities = [0.008, 0.009, 0.01, 0.011, 0.012,  0.013, 0.014, 0.015, 0.016, 0.017, 0.018]
m.run_tests(3000, 1, 7, probabilities)
    #k=8
probabilities = [0.008, 0.0009, 0.01, 0.011, 0.0115, 0.012, 0.0125, 0.013, 0.0135, 0.014, 0.015, 0.016, 0.0165, 0.017, 0.0175, 0.018, 0.019]
m.run_tests(3000, 1, 8, probabilities)
