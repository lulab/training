# define a class to hold each read
class Fastq(object):

    def __init__(self, name, seq, qual):
        self.name = name
        self.seq = seq
        self.qual = qual

    def trimDetail(self, separator=" "):
        self.temp = self.name.split(separator)
        del(self.temp[-1])
        return separator.join(self.temp)


def readFastq(infile):
    with open(infile,mode="r") as f:
        while True:
            name = f.readline().strip()
            if not name:
                break
            seq = f.readline().strip()
            name2 = f.readline().strip()
            qual = f.readline().strip()
            yield Fastq(name, seq, qual)

if __name__ == "__main__":

    in1 = '4OH1_5M_a_sub.fastq'
    in2 = '4OH1_5M_b_sub.fastq'
    seq1_dict = {}
    seq2_dict = {}
    seq1 = readFastq(in1)
    seq2 = readFastq(in2)
    s1_finished = False
    s2_finished = False

    with open(in1, "w") as reads1:
        fieldnames = ['SequenceID', 'SequenceA', 'SequenceB','QualityA','QualityB']
        writer = csv.DictWriter(reads1, fieldnames=fieldnames)
        writer.writeheader()
        while not (s1_finished and s2_finished):
            # iterate through reads (rows)
            try:
                s1 = seq1.next()
            except:
                s1_finished = True
            try:
                s2 = seq2.next()
            except:
                s2_finished = True
            # add next read to fastq1 dictionary and fastq2 dictionary
            if not s1_finished:
                seq1_dict[s1.trimDetail()] = s1

            if not s2_finished:
                seq2_dict[s2.trimDetail()] = s2

            # compare and write out reads in both dictionaries
            if not s1_finished and s1.trimDetail() in seq2_dict:
                writer.writerow({'SequenceID':seq1_dict[s1.trimDetail()].name,
                    'SequenceA':seq1_dict[s1.trimDetail()].seq,
                    'SequenceB':seq2_dict[s1.trimDetail()].seq,
                    'QualityA':seq1_dict[s1.trimDetail()].qual,
                    'QualityB':seq2_dict[s1.trimDetail()].qual})
                seq1_dict.pop(s1.trimDetail())
                seq2_dict.pop(s1.trimDetail())

            if not s2_finished and s2.trimDetail() in seq1_dict:
                writer.writerow({'SequenceID':seq1_dict[s2.trimDetail()].name,
                    'SequenceA':seq1_dict[s2.trimDetail()].seq,
                    'SequenceB':seq2_dict[s2.trimDetail()].seq,
                    'QualityA':seq1_dict[s2.trimDetail()].qual,
                    'QualityB':seq2_dict[s2.trimDetail()].qual})
                seq1_dict.pop(s2.trimDetail())
                seq2_dict.pop(s2.trimDetail())