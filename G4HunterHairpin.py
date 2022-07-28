#/bin/python
########################################################################
"""
    <G4Hunter - a program to search quadruplex-forming regions in DNA.>
    Copyright (C) <2012>  <Bedrat amina  supervised by Dr.Jean Louis Mergny>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
########################################################################
import os, re, sys, getopt
import time
import shutil
import numpy as np
from matplotlib import pyplot
import matplotlib.pyplot as plt
from Bio import SeqIO


def main(argv):
    if not argv:
        sys.stdout.write("Sorry: you must specify at least an argument\n")
        sys.stdout.write("More help avalaible with -h or --help option\n")
        sys.exit(1)

    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:w:s:", ["help", "ifile=", "ofile="])
    except getopt.GetoptError:
        print '\033[1m' + 'python G4Hunter.py -i <inputfile> -o <outputrepository> -w <window> -s <score threshold>\n' + '\033[0;0m'
        sys.exit(1)

    for opt, arg in opts:
        if opt in ('-h', "--help"):
            print  '\033[1m' + '\t ----------------------' + '\033[0;0m'
            print  '\033[1m' + '\n\t  Welcome To G4Hunter :' + '\033[0;0m'
            print  '\033[1m' + '\t ----------------------' + '\033[0;0m'

            print '\n G4Hunter takes into account G-richness and G-skewness of a given sequence and gives a quadruplex propensity score as output.'
            print 'To run G4Hunter use the commande line: \n'
            print  '\033[1m' + 'python G4Hunter.py -i <inputfile> -o <outputrepository> -w <window> -s <score threshold>\n' + '\033[0;0m'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-w", "--wind"):
            window = arg
        elif opt in ("-s", "--score"):
            score = arg

    return inputfile, outputfile, int(window), float(score)


# calcule le score de chaque base dans un fichier qui peut contenir plusieur sequences fasta
# le fichier doit comporter la nom de la seq + la sequence en une seul ligne donc pas de \n ou \r
class Soft(object):

    def __init__(self):
        pass

    def ReadFile(self, Filein):
        ListSeq, LHeader = [], []
        for record in SeqIO.parse(Filein, "fasta"):
            LHeader.append(record.id)
            ListSeq.append(record.seq)
        return LHeader, ListSeq

    def ReadSeq(self, Seqs):
        ListSeq, LHeader = [], []
        for record in SeqIO.parse(Seqs, "fasta"):
            LHeader.append(record.id)
            ListSeq.append(record.seq)
        return LHeader, ListSeq

    def GFinder(self, Filein, k):
        LHeader, ListSeq = self.ReadFile(Filein)
        LSeq, LMark, LNumber, LScoreSeq, SeqLine = [], [], [], [], ""
        for i in range(len(ListSeq)):
            liste, utter_left_position, utter_right_position, mark, Sequence = self.G_scoring(ListSeq[i])
            LSeq.append(Sequence)
            LMark.append(mark)
            LScoreSeq.append(self.CalScore(liste, k))
            LNumber.append(liste)
        return LScoreSeq, LSeq, LMark, LNumber, LHeader


    def BP(self, line): #rekurzivno iskanje hairpinov, primerjava med razlicnimi hairpini in iskanje najboljsega
        base_pairs = [("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("A", "T"), ("T", "A"), ("U", "G"), ("G", "U"), ("G","T"), ("T, G")] #GT bazni par zaradi GU pri RNA zaporedjih v NCBI, za DNA ga daj ven
        #in se hoogsteen
        if len(line) == 0:
            return ("", [], 0)
        #elif len(line) == 1:
            #return (".", [0.2], 0.2)
        #elif len(line) == 2:
            #return ("..", [0.2, 0.2], 0.4)
        elif len(line) == 3:
            return ("...", [0.2, 0.2, 0.2], 0.6)
        # elif len(line)=2 and (x,z) in base_pairs:
        elif len(line) == 4:
            return ("....", [0.2, 0.2, 0.2, 0.2], 0.8)

        else:
            x = line[0]
            y = line[1:-1]
            z = line[-1]

            if (x, z) in base_pairs:
                #if len(y) == 0:
                    #return ("()", [-4, -4], -8)  # negativen score da nikoli ni te situacije
                result = self.BP(y)
                return ("(" + result[0] + ")", [0.6] + result[1] + [0.6], 1.2 + result[2])
            else:
                result1 = self.BP(x + y)
                result2 = self.BP(y + z)
                if result1[2] > result2[2]:
                    return (result1[0] + ".", result1[1] + [-1.2], result1[2] -1.2)
                else:
                    # result=result2
                    return ("." + result2[0], [1.2] + result2[-2], result2[2] - 1.2)

    def G_scoring(self, seq): #rekurzivna reformulacija G4Hunter
        if len(seq) == 1:
            if seq == "G":
                score_list = [1]
                counter = 1
                # score_total=1
                utter_left_position = 0
                utter_right_position = 0
                mark = "G"

            elif seq == "C":
                score_list = [-1]
                counter = 0
                # score_total=1
                utter_left_position = -1
                utter_right_position = -1
                mark = "C"

            #elif seq == "A" or seq == "U" or seq == "T":
            else:
                score_list = [0]
                # score_total=1
                counter = 0
                utter_left_position = -1
                utter_right_position = -1
                mark = "."

            return score_list, utter_left_position, utter_right_position, mark, seq
        #razpolavljamo zaporedje do dna rekurzije
        middle = len(seq) // 2
        # if middle
        seq1 = seq[:middle]
        seq2 = seq[middle:]
        score_list1, utter_left_position1, utter_right_position1, mark1, seq1 = self.G_scoring(seq1)
        score_list2, utter_left_position2, utter_right_position2, mark2, seq2 = self.G_scoring(seq2)

        # preverimo ali je v seq1 samo desni G otok
        if utter_left_position1 >= 0:
            utter_left_position = utter_left_position1
        elif utter_right_position2 >= 0:
            utter_left_position = len(score_list1) + utter_left_position2
        else:
            utter_left_position = utter_left_position2

        if utter_right_position2 >= 0:
            utter_right_position = len(score_list1) + utter_right_position2
        else:
            utter_right_position = utter_right_position1

        utter_right_base1 = seq1[len(seq1) - 1]
        utter_left_base2 = seq2[0]
        utter_right_score1 = score_list1[len(seq1) - 1]
        utter_left_score2 = score_list2[0]
        new_island_score = utter_right_score1 + utter_left_score2

        mark = mark1 + mark2

        i = 0
        score_list = score_list1 + score_list2

        # zdruzevanje G otokov, ko pride do stika sosednjih G-jev
        if utter_right_base1 == "G" and utter_left_base2 == "G":
            G_counter1 = 1
            G_counter2 = 1
            while utter_right_position1 - G_counter1 > 0 and seq1[utter_right_position1 - G_counter1] == "G":
                G_counter1 += 1
            while utter_left_position2 + G_counter2 < len(seq2) and seq2[utter_left_position2 + G_counter2] == "G":
                G_counter2 += 1
            G_counter = G_counter1 + G_counter2
            #da je scoring enega G drugacen sistem kot pri vec G
            for i in range(G_counter):
                score_list[utter_right_position1 - G_counter1 + 1 + i] = G_counter+2

        # zdruzevanje C otokov
        elif utter_right_base1 == "C" and utter_left_base2 == "C":

            for i in range(-new_island_score):
                score_list[len(seq1) - 1 + utter_right_score1 + 1 + i] = new_island_score
        else: #iskanje hairpinov med G - primerjava vseh tockovanj
            for window_size in range (9,15):
                if len(seq1)>window_size:
                    for window_scan in range(len(seq1)-window_size, len(seq1)):
                        if seq[window_scan-1] == "G" and seq[window_scan+window_size] == "G":
                            seq3 = seq[window_scan:window_scan + window_size]
                            hairpin = self.BP(seq3)
                            hairpin_score = hairpin[1]
                            hairpin_mark = hairpin[0]
                            sum_score = 0
                            for i in range(window_size):
                                sum_score += score_list[window_scan + i]
                                i += 1
                            #print(sum_score)

                            if hairpin[2] > sum_score and hairpin[2] > 3:
                                score_list = score_list[:window_scan] + hairpin_score + score_list[window_scan+window_size:]
                                mark = mark[:window_scan] + hairpin_mark + mark[window_scan+window_size:]

            #seq3 = seq[(utter_right_position1 + 1):(len(seq1) + utter_left_position2)]
            #if 20 > len(seq3) >= 7:
               # hairpin = self.BP(seq3)
               # hairpin_score = hairpin[1]
               # hairpin_mark = hairpin[0]
               # if hairpin[2] > 3:
               #     score_list = score_list[:utter_right_position1 + 1] + hairpin_score + score_list[len(seq1) + utter_left_position2:]
               #      mark = mark[:utter_right_position1 + 1] + hairpin_mark + mark[len(seq1) + utter_left_position2:]

        # hairpin iskanje ko 6<len<33

        return score_list, utter_left_position, utter_right_position, mark, seq

    def CalScore(self, liste, k):
        Score_Liste = []
        # calcule de la moynne des scores pour toutes les sequences - les derniers k bases
        for i in range(len(liste) - (k - 1)):
            # print (len(liste)), i, k
            j, Sum = 0, 0
            while (j < k):
                # print j, i
                Sum = Sum + liste[i]
                j = j + 1
                i = i + 1
            Mean = Sum / float(k)
            Score_Liste.append(Mean)
        return Score_Liste

    ###################################################
    ###################################################
    def plot2(self, liste, repert, i):
        # make a little extra space between the subplots
        plt.subplots_adjust(wspace=1.0)
        dt = 1
        t = np.arange(0, len(liste), dt)
        figure = plt.figure()
        plt.plot(t, liste, 'b-')
        plt.xlim(0, len(liste))
        # figure.suptitle('Score of window nucleotides', fontsize=16)
        plt.xlabel('Position (ntS)')
        plt.ylabel('Score')
        plt.grid(True)
        figure.savefig(repert + 'Score_plot_' + i + '.pdf', dpi=figure.dpi)

    """ 
    ###################################################
    ###################################################
    """

    def GetG4(self, line, mark, fileout, liste, Window, k, header, Len):
        LG4 = []
        SEQ = ">" + header + "\n Start \t End \t Dot bracket \t\t Length \t Score \t Sequence \n"
        fileout.write(SEQ)
        for i in range(len(liste)):
            if (liste[i] >= float(Window) or liste[i] <= - float(Window)):
                seq = line[i:i + k]
                dots = mark[i:i+k]
                LG4.append(i)
                self.Write(fileout, i, k, 0, 0, seq, dots, k, liste[i])
                fileout.write("\n")
        return LG4

    def WriteSeq(self, line, mark, fileout, liste, LISTE, header, F, Len):
        i, k, I = 0, 0, 0
        a = b = LISTE[i]
        MSCORE = []
        SEQ = ">" + header + "\nStart\tEnd\tDot bracket\t\tLength\tScore\tNBR\tSequence\n"
        fileout.write(SEQ)
        if (len(LISTE) > 1):
            c = LISTE[i + 1]
            while (i < len(LISTE) - 2):
                if (c == b + 1):
                    k = k + 1
                    i = i + 1
                else:
                    I = I + 1
                    seq = line[a:a + F + k]
                    dots = mark[a:a + F + k]
                    liste2, utter_left_position, utter_right_position, mark2, sequence = self.G_scoring(seq)
                    self.Write(fileout, a, k, F, 0, seq, dots, len(seq), round(np.mean(liste2), 2))
                    MSCORE.append(abs(round(np.mean(liste2), 2)))
                    fileout.write("\n")
                    k = 0
                    i = i + 1
                    a = LISTE[i]
                b = LISTE[i]
                c = LISTE[i + 1]
            I = I + 1
            seq = line[a:a + F + k + 1]
            dots = mark[a:a + F + k +1]
            liste2, utter_left_position, utter_right_position, mark2, sequence = self.G_scoring(seq)
            self.Write(fileout, a, k, F, 1, seq, dots, len(seq), round(np.mean(liste2), 2))
            MSCORE.append(abs(round(np.mean(liste2), 2)))
            fileout.write("\t")
            fileout.write(str(I))
            fileout.write("\n")
            # dans le cas ou il ya une seul sequence donc pas de chevauchement
        else:
            I = I + 1
            seq = line[a:a + F]
            dots = mark[a:a + F]
            self.Write(fileout, a, 0, F, 0, seq, dots, len(seq), liste[a])
            MSCORE.append(abs(liste[a]))
            fileout.write("\t")
            fileout.write(str(I))
            fileout.write("\n")
        return MSCORE

    def Write(self, fileout, i, k, F, X, seq, mark, long, score):
        LINE = str(i) + " \t " + str(i + k + F + X) + " \t " + str(mark) + " \t\t" + str(long) + " \t " + str(score) + " \t " + str(seq)
        fileout.write(LINE)

    # Len dans le cas ou la sequence fini avec des ----- donc il yaura une erreur


if __name__ == "__main__":
    try:
        inputfile, outputfile, window, score = main(sys.argv[1:])
        fname = inputfile.split("/")[-1]
        name = fname.split(".")

    except ValueError:
        print '\033[1m' + "\n \t Oops! invalide parameters  \n" + '\033[0;0m'
        print "--------------------------------------------------------------------\n"
        sys.exit()
    except UnboundLocalError:
        print '\033[1m' + "\n \t Oops! invalide parameters  \n" + '\033[0;0m'
        print "--------------------------------------------------------------------\n"
        sys.exit()

    OPF = os.listdir(outputfile)
    flag = False
    for dir in OPF:
        DIR = "Results_" + str(name[0])
        if dir == DIR:
            print "true", DIR
            flag = True
    if flag == True:
        shutil.rmtree(outputfile + "/" + DIR + "/")
        os.makedirs(outputfile + "/" + DIR + "/", mode=0777)  #
        print '\033[1m' + "\n \t Re-evaluation of G-quadruplex propensity with G4Hunter " + '\033[0;0m'
        print "\n#####################################"
        print "#    New Results directory Created  #"
        print "#####################################\n"
    else:
        os.makedirs(outputfile + "/" + DIR + "/", mode=0777)  #
        print "\n########################################################################"
        print "#                            Results directory Created                 #"
        print "########################################################################\n"

    # ================================================================
    plot = []
    fname = inputfile.split("/")[-1]
    filefasta = fname.split(".")
    filein = open(inputfile, "r")
    print "\n Input file:", '\033[1m' + filefasta[0] + '\033[0;0m'
    # repertoire des fichiers de sortie

    Res1file = open(outputfile + "/" + DIR + "/" + filefasta[0] + "-W" + str(window) + "-S" + str(score) + ".txt", "w")
    Res2file = open(outputfile + "/" + DIR + "/" + filefasta[0] + "-Merged.txt", "w")
    # =========================================

    startTime = time.time()

    soft1 = Soft()
    ScoreListe, DNASeq, DNAMark, NumListe, HeaderListe = soft1.GFinder(filein, window)
    for i in range(len(DNASeq)):
        G4Seq = soft1.GetG4(DNASeq[i], DNAMark[i], Res1file, ScoreListe[i], float(score), int(window), HeaderListe[i],
                            len(NumListe[i]))
        if (len(G4Seq) > 0):
            MSCORE = soft1.WriteSeq(DNASeq[i], DNAMark[i], Res2file, ScoreListe[i], G4Seq, HeaderListe[i], int(window),
                                    len(NumListe[i]))
    # plot.append(MSCORE)
    """
    malist, alllist=[], []
    #print ScoreListe[0]
    for jj in range (len(ScoreListe[0])):
        cc, mean=0, 0
        for kk in range(len(ScoreListe)-1):
            #print kk, jj, ScoreListe[kk][jj]
            cc+=ScoreListe[kk][jj]
        mean=cc/len(ScoreListe)
        alllist.append(mean)
        if abs(mean) >=score :
            malist.append(mean)
        else:
            malist.append(0)

    #soft1.plot2(ScoreListe[0], outputfile +"/Results/")
    soft1.plot2(malist, outputfile +"/"+DIR+"/", "sc")
    soft1.plot2(alllist, outputfile +"/"+DIR+"/", "all")
    """
    filein.close()
    fin = time.time()

    print "\n Results files and Score Figure are created in:   "  # ,fin-startTime, "secondes"
    print '\033[1m' + outputfile, "/", DIR, "/", "\n " + '\033[0;0m'

    Res1file.close()
    Res2file.close()


