from django.urls import path, include

from scipy import signal

from .detect_peaks import detect_peaks
import matplotlib.pyplot as plt
import numpy as np
from .findMax import findMax
import pickle


def classify(data):
    # dataFile = open("data2.csv", "w")
    count = 0
    for num in range(1):

        # if num == 365:
        #     continue
        # filename = "Nsignal%d.csv" % (num + 1)
        # file = open(filename, "r")
        #
        # myFile = file.read()
        new1 = data.replace('[', '')
        new2 = new1.replace(']', '')
        new3 = new2.split('\n')
        new = new3[0].split(",")
        ecg = []
        for i in range(len(new) - 1):
            ecg.append(float(new[i]))

        # print(ecg[:])
        # plt.figure(1)
        # plt.plot(ecg)

        qrs_c = []
        qrs_i = []
        SIG_LEV = 0
        nois_c = []
        nois_i = []
        delay = 0
        skip = 0
        not_nois = 0
        selected_RR = []
        m_selected_RR = 0
        mean_RR = 0
        qrs_i_raw = []
        qrs_amp_raw = []
        ser_back = 0
        test_m = 0
        SIGL_buf = []
        NOISL_buf = []
        THRS_buf = []
        SIGL_buf1 = []
        NOISL_buf1 = []
        THRS_buf1 = []
        fs = 200
        ecg[:] = [x - (sum(ecg) / len(ecg)) for x in ecg]
        # ecg_n = np.divide(ecg, max(ecg))
        # print(ecg)
        Wn = 12 * 2 / fs
        N = 3
        [a, b] = signal.butter(3, Wn, 'low')
        # print(a)
        # print(b)
        ecg_1 = signal.filtfilt(a, b, ecg)
        ecg_1 = np.divide(ecg_1, max(abs(ecg_1)))

        # plt.figure(2)
        # plt.plot(ecg_1)


        Wn = 5 * 2 / fs
        N = 3
        [a, b] = signal.butter(3, Wn, 'high')
        ecg_h = signal.filtfilt(a, b, ecg_1)
        ecg_h = np.divide(ecg_h, max(abs(ecg_h)))

        # plt.figure(3)
        # plt.plot(ecg_h)

        o = [1, 2, 0, -2, -1]
        b = [x * (fs / 8) for x in o]

        ecg_d = signal.filtfilt(b, 1, ecg_h)
        ecg_d = np.divide(ecg_d, max(abs(ecg_d)))

        # plt.figure(4)
        # plt.plot(ecg_d)

        ecg_s = np.square(ecg_d)

        # plt.figure(5)
        # plt.plot(ecg_s)

        window = np.round(fs * 0.15)
        weights = np.repeat(1.0, window) / window
        ecg_m = np.convolve(ecg_s, weights, 'valid')
        ecg_m = np.divide(ecg_m, max(abs(ecg_m)))

        # plt.figure(6)
        locs = detect_peaks(ecg_m, mpd=40, show=False)
        # locs=signal.find_peaks_cwt(ecg_m, np.arange(1.0, 40))
        pks = np.arange(1.0, len(locs) + 1)
        for l in range(len(locs)):
            pks[l] = ecg_m[locs[l]]

        # plt.plot(ecg_m)

        THR_SIG = max(ecg_m[1:2 * fs]) * 1 / 3
        THR_NOISE = np.mean(ecg_m[1:2 * fs]) * 1 / 2
        SIG_LEV = THR_SIG
        NOISE_LEV = THR_NOISE

        THR_SIG1 = max(ecg_h[1:2 * fs]) * 1 / 3
        THR_NOISE1 = np.mean(ecg_h[1:2 * fs]) * 1 / 2
        SIG_LEV1 = THR_SIG1
        NOISE_LEV1 = THR_NOISE1

        for i in range(len(locs)):
            if (locs[i] - np.round(0.15 * fs) >= 1) and (locs[i] <= len(ecg_h)):
                x_i, y_i = findMax((ecg_h[locs[i]:locs[i] + round(0.15 * fs)]))
                # print(ecg_h[locs[i]-round(0.15*fs):locs[i]])

            else:
                if i == 0:
                    x_i, y_i = findMax(ecg_h[:locs[i]])
                    ser_back = 1

                elif locs[i] >= len(ecg_h):
                    x_i, y_i = findMax(ecg_h[locs[i] - round(0.15 * fs):])
            if len(qrs_c) >= 9:
                diffRR = np.diff(qrs_i[len(qrs_i) - 8:])
                mean_RR = np.mean(diffRR)
                comp = qrs_i[len(qrs_i) - 1] - qrs_i[len(qrs_i) - 1]

                if comp <= 0.92 * mean_RR or comp >= 1.16 * mean_RR:
                    THR_SIG = 0.5 * THR_SIG
                    THR_SIG1 = 0.5 * THR_SIG1

                else:
                    m_selected_RR = mean_RR

            if m_selected_RR:
                test_m = m_selected_RR
            elif mean_RR and m_selected_RR == 0:
                test_m = mean_RR
            else:
                test_m = 0

            if test_m:
                if (locs[i] - qrs_i[len(qrs_i) - 1]) >= round(1.66 * test_m):
                    locs_temp, pks_temp = findMax(
                        ecg_m[qrs_i[len(qrs_i) - 1] + round(0.200 * fs):locs[i] - round(0.200 * fs)])

                    locs_temp = qrs_i[len(qrs_i) - 1] + round(0.200 * fs) + locs_temp - 1

                    if pks_temp > THR_NOISE:
                        qrs_c.append(pks_temp)
                        qrs_i.append(locs_temp)

                        if locs_temp <= len(ecg_h):
                            x_i_t, y_i_t = findMax(ecg_h[locs_temp - round(0.150 * fs):locs_temp])

                        else:
                            x_i_t, y_i_t = findMax(ecg_h[locs_temp - round(0.150 * fs):])

                        if y_i_t > THR_NOISE1:
                            qrs_i_raw.append(locs_temp - round(0.150 * fs) + (x_i_t - 1))
                            qrs_amp_raw.append(y_i_t)
                            SIG_LEV1 = 0.25 * y_i_t + 0.75 * SIG_LEV1

                        not_nois = 1
                        SIG_LEV = 0.25 * pks_temp + 0.75 * SIG_LEV
                else:
                    not_nois = 0

            if pks[i] >= THR_SIG:
                if len(qrs_c) >= 3:
                    if locs[i] - qrs_i[len(qrs_i) - 1] <= round(0.36 * fs):
                        # print(qrs_i[len(qrs_i)-1-round(0.075*fs):qrs_i[i]])

                        slope1 = np.mean(np.diff(ecg_m[locs[i] - round(0.075 * fs):locs[i]]))
                        slope2 = np.mean(np.diff(ecg_m[qrs_i[len(qrs_i) - 1] - round(0.075 * fs):qrs_i[len(qrs_i) - 1]]))

                        if abs(slope1) <= abs(0.5 * slope2):
                            nois_c.append(pks[i])
                            nois_i.append(locs[i])
                            skip = 1
                            NOISE_LEV1 = 0.125 * y_i + 0.875 * NOISE_LEV1
                            NOISE_LEV = 0.125 * pks[i] + 0.875 * NOISE_LEV
                        else:
                            skip = 0
                if skip == 0:
                    qrs_c.append(pks[i])
                    qrs_i.append(locs[i])

                    # print(y_i)
                    # print(THR_SIG1)
                    if y_i >= THR_SIG1:
                        if ser_back:
                            qrs_i_raw.append(x_i)
                        else:
                            qrs_i_raw.append(locs[i] - round(0.150 * fs) + (x_i - 1))
                    qrs_amp_raw.append(y_i)
                    SIG_LEV1 = 0.125 * y_i + 0.875 * SIG_LEV1
                SIG_LEV = 0.125 * pks[i] + 0.875 * SIG_LEV

            elif THR_NOISE <= pks[i] < THR_SIG:
                NOISE_LEV1 = 0.125 * y_i + 0.875 * NOISE_LEV1
                NOISE_LEV = 0.125 * pks[i] + 0.875 * NOISE_LEV

            elif pks[i] < THR_NOISE:
                nois_c.append(pks[i])
                nois_i.append(locs[i])

                NOISE_LEV1 = 0.125 * y_i + 0.875 * NOISE_LEV1
                NOISE_LEV = 0.125 * pks[i] + 0.875 * NOISE_LEV

            if NOISE_LEV != 0 or SIG_LEV != 0:
                THR_SIG = NOISE_LEV + 0.25 * (abs(SIG_LEV - NOISE_LEV))
                THR_NOISE = 0.5 * THR_SIG

            if NOISE_LEV1 != 0 or SIG_LEV1 != 0:
                THR_SIG1 = NOISE_LEV1 + 0.25 * (abs(SIG_LEV1 - NOISE_LEV1))
                THR_NOISE1 = 0.5 * THR_SIG1

            SIGL_buf.append(SIG_LEV)
            NOISL_buf.append(NOISE_LEV)
            THRS_buf.append(THR_SIG)

            SIGL_buf1.append(SIG_LEV1)
            NOISL_buf1.append(NOISE_LEV1)
            THRS_buf1.append(THR_SIG1)

            skip = 0
            not_nois = 0
            ser_back = 0

        for asd in range(len(qrs_i_raw)):
            qrs_i_raw[asd] += 31

        new = np.arange(1.0, len(qrs_i_raw) + 1)
        for zxc in range(len(qrs_i_raw)):
            new[zxc] = ecg_1[qrs_i_raw[zxc]]

        RRInterval = []
        for value in range(len(qrs_i_raw) - 2):
            RRInterval.append(qrs_i_raw[value + 1] - qrs_i_raw[value])

        # print("hello")

        # for asd in range(len(qrs_i_raw)):
        #    qrs_i_raw[asd] += 31

        new = np.arange(1.0, len(qrs_i_raw) + 1)
        for zxc in range(len(qrs_i_raw)):
            new[zxc] = ecg_h[qrs_i_raw[zxc]]

        RRInterval = []
        for value in range(len(qrs_i_raw) - 2):
            RRInterval.append(qrs_i_raw[value + 1] - qrs_i_raw[value])

        avgRPeak = np.average(new)
        stdRPeak = np.std(new)
        avgRRInterval = np.average(RRInterval)
        stdRRInterval = np.std(RRInterval)

        power = (sum(np.square(ecg)) * (1 / len(ecg)))

        print("Features: %f, %f, %f, %f, %f" % (avgRPeak, avgRRInterval, stdRPeak, stdRRInterval, power))

        randomForest = pickle.load(open("modelrand2", 'rb'))
        predictionrand = randomForest.predict([[avgRPeak, stdRPeak, avgRRInterval, stdRRInterval, power]])
        # if num <=71:
        #    if prediction == 1:
        #        count=count+1
        # else:
        if predictionrand == 1:
            count = count + 1

    return str(predictionrand), avgRPeak
#    print(str(num + 1) + ": " + str(predictionrand) + "\n")
"""
    plt.figure()
    plt.plot(ecg[500:1700])
    plt.figure(2)
    plt.plot(ecg_1[500:1700])
    plt.figure(3)
    plt.plot(ecg_h[500:1700])
    plt.figure(4)
    plt.plot(ecg_d[500:1700])
    plt.figure(5)
    plt.plot(ecg_s[500:1700])
    plt.figure(6)
    plt.plot(ecg_m)
    plt.figure(7)
    plt.figure()
    ecg_n = np.divide(ecg, max(ecg[:]))
    plt.plot(ecg_n)
    plt.scatter(qrs_i_raw, new, color="r")
plt.show()
print(count)
"""
