# This is a sptxt file parser
# read the file and return the spectrum.



class LibItem:
    def __init__(self):
        self.seq = ''
        self.charge = ''
        self.Name = ''
        self.LibID = ''
        self.Status = ''
        self.PrecursorMZ = ''
        self.FullName = ''
        self.MW=''
        self.Comment = ''
        self.NumPeaks=''
        self.mz=[]
        self.intensity=[]
        self.annotations=[]
    def print(self):
        print('seq={},charge={},name={},libID={},status={},precursormz={},'
              'fullname={},MW={},Comment={}'.format(self.seq, self.charge,self.Name, self.LibID,
                                                    self.Status, self.PrecursorMZ, self.FullName, self.MW,self.Comment),end='')
        print('\npeaks:\nmz={}\nintensity={}\nannotations={}'.format(self.mz, self.intensity,self.annotations))


class SpTXT:
    def __init__(self):
        self.speclib=[]
        self.seqchargetoidx={}
        self.precursormzs=[]


    def parse(self, filename):
        self.sptxtfilename = filename
        fid = open(filename, 'r')
        lines = fid.readlines()
        fid.close()

        self.speclib = []
        i=0
        while i <len(lines):
            # print(lines[i],end='')
            l = lines[i]
            i+=1
            if len(l)>0:
                if l[0]=='#':
                   print('ignore comment',end='\n')
                if  'Name:' in l and l[:5]=='Name:':
                    if(len(self.speclib)%1000==0): print('._.{}'.format(len(self.speclib)),end='',flush=True)
                    sequence_charge = l.split()[1]
                    self.speclib.append(LibItem())
                    self.speclib[-1].Name = sequence_charge
                    [seq, charge] = sequence_charge.split('/')
                    self.speclib[-1].charge = charge
                    self.speclib[-1].seq = seq

                    # speclib[-1].print()

                if 'LibID:' in l and l[:6] == 'LibID:':
                    libid = l.split()[1]
                    self.speclib[-1].LibID = libid

                if 'MW:' in l and l[:3] == 'MW:':
                    MW = l.split()[1]
                    self.speclib[-1].MW = MW

                if 'PrecursorMZ:' in l and l[:12] == 'PrecursorMZ:':
                    PrecursorMZ = l.split()[1]
                    self.speclib[-1].PrecursorMZ = PrecursorMZ

                if 'FullName:' in l and l[:9] == 'FullName:':
                    FullName = l.split()[1]
                    self.speclib[-1].FullName = FullName



                if 'NumPeaks:' in l and l[:9] == 'NumPeaks:':
                    numpeaks = l.split()[1]
                    self.speclib[-1].NumPeaks = numpeaks
                    self.speclib[-1].mz=[0.0]*int(numpeaks)
                    self.speclib[-1].intensity=[0.0]*int(numpeaks)
                    self.speclib[-1].annotations=['']*int(numpeaks)


                    # get peaks
                    for j in range(int(numpeaks)):
                        pk = lines[i+j]
                        # print(pk.split())
                        mz, inten, annotation, freq, info = pk.split()
                        self.speclib[-1].mz[j]=float(mz)
                        self.speclib[-1].intensity[j] = float(inten)
                        self.speclib[-1].annotations[j] = annotation

                # if 'LibID:' in l and l[:6] == 'LibID:':
                #     libid = l.split()[1]
                #     speclib[-1].LibID = libid
        # end of while
        print('._.{}'.format(len(self.speclib)), end='\n', flush=True)

        for k in range(len(self.speclib)):
            item = self.speclib[k]
            self.seqchargetoidx[item.Name]=k

        print('Index built for quick access!')

    def getLibItemIdx(self,seqcharge):
        if seqcharge in self.seqchargetoidx:
            return self.seqchargetoidx[seqcharge]
        else:
            print('Seq charge \"{}\"not found'.format(seqcharge))
            return -1

    def getLibItem(self, seqcharge):
        idx = self.getLibItemIdx(seqcharge)
        return self.getLibItemOnIDX(idx)

    def getLibItemOnIDX(self, idx):
        if idx == -1 or idx >= len(self.speclib):
            return LibItem()
        else:
            return self.speclib[idx]

    def counts(self):
        return len(self.speclib)

    def getIdx(self, precursorMZ, tolerance):
        if len(self.precursormzs) == 0:
            for i in range(len(self.speclib)):
                self.precursormzs.append(float(self.speclib[i].PrecursorMZ))

        l=[]
        for i in range(len(self.precursormzs)):
            if abs(precursorMZ-self.precursormzs[i])<tolerance:
                l.append(i)
        return l
