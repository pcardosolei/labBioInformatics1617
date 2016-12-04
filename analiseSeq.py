"""
  Parte BioPython
"""

from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio import ExPASy
from Bio import SwissProt


"""
  Imports genericos
"""
import os
import re
import time

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

  def addProteinID(self,ident):
    self.id = ident

  def addProduct(self,prod):
    self.product = prod
""" ###########

  Aspetos gerais na analise da sequencia e features presentes na NCBI
""" ###########

"""
  Ler a sequencia do NCBI

  parametros - utilizador , posicao inicial da sequencia , posicao final da sequencia
"""
def readNCBI(user,startSeq,stopSeq):
  Entrez.email = user
  handle = Entrez.efetch(db="nucleotide",id="NC_002942.5",rettype="gbwithparts",retmode="text",seq_start=startSeq,seq_stop=stopSeq)
  return handle

# escrever a sequencia para um ficheiro

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
        if 'protein_id' in qualifiers.keys():
          aux.addProteinID(qualifiers['protein_id'][0])
        if 'product' in qualifiers:
          aux.addProduct(qualifiers['product'][0])
        tag = qualifiers['db_xref'][0]
        dic[tag] = aux
    else:
      pass
  return dic

"""
  Funcao para obter valores adicionais de cada um dos genes -- ja nao me lembro onde xD
"""
def getInfoFeature():
  pass


def getProteinFeatures():
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
  for tag in dic:
    if hasattr(dic[tag], 'seq'):
      ofile = open("fastas/"+tag+".fasta", "w")
      ofile.write(">" + tag + "\n" +dic[tag].seq + "\n")
      #do not forget to close it
      ofile.close()
    else:
      pass


#http://biopython.org/DIST/docs/api/Bio.Blast.NCBIWWW-module.html
# ESTUDAR ESTE LINK PARA VER OS PARAMETROS QUE INTERESSAM

def fastaToBlastFile(fasta):
  ofile = open(fasta,"r")
  seq = ""
  for line in ofile:
    seq += line
  result_handle = NCBIWWW.qblast("blastp", "nr",seq)
  name = os.path.basename(fasta)
  blastFile = open("blasts/"+os.path.splitext(name)[0]+".xml","w")
  blastFile.write(result_handle.read())
  blastFile.close()
  result_handle.close()
  ofile.close()




def blastAllSeq(folder):
  sum = 0
  fastas = []
  fastas += [each for each in os.listdir(folder) if each.endswith('.fasta')]

  for fasta in fastas:
    fastaToBlastFile(fasta)
    sum += 1
    time.sleep(300) #pausar a funcao por 5 minutos para nao estourar
  print("Foram tratadas {} sequencias".format(sum))

"""
  SWISS PROT Function
"""

# not working não apanha a tag da legionella
def getSwissProtInfo(tag):
  handle = ExPASy.get_sprot_raw(tag)
  html_results = handle.read()
  ofile = open("teste.html","wb")
  ofile.write(html_results)
  ofile.close()


"""
  Create folders
"""
def createFolders(directories):
  for directory in directories:
    if not os.path.exists(directory):
      os.makedirs(directory)

"""
  Guardar informacao genes
"""

def saveGenes(dictionary):
  ofile = open("genes.txt","w")
  for gene in dictionary:
    txt = stringInfo(dictionary[gene])
    ofile.write(gene+"\n")
    ofile.write(txt)
  ofile.close()

def stringInfo(feature):
  text = ""
  if hasattr(feature,"locus_tag"):
    text += "Locus_tag: "+ feature.locus_tag + "\n"
  if hasattr(feature,"function"):
    text += "Function: "+ ''.join(feature.function) + "\n"
  if hasattr(feature,"note"):
    text += "Note: " + feature.note + "\n"
  if hasattr(feature,"id"):
    text += "Protein_id: " + feature.id + "\n"
  if hasattr(feature,"product"):
    text += "Product: " + feature.id + "\n"
  text += " \n"
  return text

"""
  Testes que ainda necessitam de ser vistos
"""
def blastProtein(seq):
  result_handle = NCBIWWW.qblast("blastp","nt",seq)
  result_handle = open("my_blast.xml")



def readBlastProtein():
  result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
  print(result_handle) #isto dá como resultado um cStringIO . tenho de por isto a ler com um IO qualquer


"""
  MAIN
"""
def main():
  #folders = ["fastas","blasts"]
  #createFolders(folders)

  #readNCBI("history.a61043@alunos.uminho.pt","2873801","3124550")
  dic = locateFeatures()
  #proteinSequencesToFasta(dic)

  saveGenes(dic)
  #fastaToBlastFile("fastas/GeneID:19834107.fasta")
  #blastAllSeq("fastas")
  #getSwissProtInfo("UP000000609")

  #readBlastProtein()
  #getUniProt()

main()
