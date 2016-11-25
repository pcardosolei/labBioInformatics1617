"""
  Parte BioPython
"""

from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio import ExPASy
from Bio import SwissProt


"""
  Imports genéricos
"""
import os

class Feature:
  def __init__(self,ref):
    self.ref = ref

  #adicionar a funcao na classe feature
  def addFunction(self,function):
    self.function = function

  def addNote(self,note):
    self.note = note

  #adicionar a locusTag na classe feature
  def addLocusTag(self,locus_tag):
    self.locus_tag = locus_tag

  #adicionar a sequencia na classe feature
  def addSeq(self,seq):
    self.seq = seq

  #e preciso testar as seguintes funcoes
  def getFunction(self):
    return self.function

  def getNote(self):
    return self.note

  def getLocusTag(self):
    return self.locus_tag

  def getSeq(self):
    return self.seq

""" ###########

  Aspetos gerais na analise da sequencia e features presentes na NCBI
""" ###########

"""
  Ler a sequencia do NCBI
"""
def readNCBI(user,startSeq,stopSeq):
  Entrez.email = user
  handle = Entrez.efetch(db="nucleotide",id="NC_002942.5",rettype="gbwithparts",retmode="text",seq_start=startSeq,seq_stop=stopSeq)
  return handle

# escrever a sequencia para um ficherio
def writeSequence(user,startSeq,stopSeq):
  handle = readNCBI(user,startSeq,stopSeq)
  record = SeqIO.read(handle,"genbank")
  SeqIO.write(records,"example.gb","genbank")

"""
  Funcao trabalha com a sequencia.
  Le as features e guarda as informacoes que sao consideradas importantes - db_ref , function e locus_tag

  Falta talvez adicionar o ID da proteina e mais cenas
"""

def locateFeatures():
  records = SeqIO.read(open("example.gb"),"genbank")
  dic = {}
  for feature in records.features:
    #vale a pena considerar a source e o STS presente na nossa proteina
    if ("CDS" or "gene") in feature.type:
      qualifiers = feature.qualifiers
      if "db_xref" in qualifiers:
        aux = Feature(qualifiers['db_xref'][0])
        if "locus_tag" in qualifiers.keys():
          aux.addLocusTag(qualifiers['locus_tag'][0])
        if "function" in qualifiers.keys():
          aux.addFunction(qualifiers['function'])
        if 'translation' in qualifiers.keys():
          aux.addSeq(qualifiers['translation'][0])
        if 'note' in qualifiers.keys():
          aux.addNote(qualifiers['note'][0])
        tag = qualifiers['db_xref'][0]
        dic[tag] = aux
    else:
      pass
  return dic

"""
  Funcao para obter valores adicionais de cada um dos genes
"""

def getInfoFeature():
  pass

"""
  UNIPROT Functions
"""
def getUniProt():
  handle = ExPASy.get_sprot_raw("P37032")
  record = SwissProt.read(handle)
  for ref in record.references:
    print(ref.title)


"""

  BLAST Functions

"""

def proteinSequencesToFasta(dic):
  sum = 0
  for tag in dic:
    if hasattr(dic[tag], 'seq'):
      ofile = open("fastas/"+tag+".fasta", "w")
      ofile.write(">" + tag + "\n" +dic[tag].seq + "\n")
      #do not forget to close it
      ofile.close()
    else:
      pass



"""
  Criar pastas
"""

def createFolders(directory):
  if not os.path.exists(directory):
    os.makedirs(directory)

"""
  This is a test
"""
def blastProtein(seq):
  result_handle = NCBIWWW.qblast("blastp","nt",seq)
  result_handle = open("my_blast.xml")



def readBlastProtein():
  result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
  print(result_handle) #isto dá como resultado um cStringIO . tenho de por isto a ler com um IO qualquer



def test():
  createFolders("fastas")

  #readNCBI("history.a61043@alunos.uminho.pt","2873801","3124550")
  dic = locateFeatures()
  proteinSequencesToFasta(dic)
  #readBlastProtein()
  #getUniProt()

test()
