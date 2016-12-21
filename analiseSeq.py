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
  def addDbRef(self,db_xref):
    self.db_xref = db_xref

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
  Ler a sequencia do NCBI da legionella

  parametros - utilizador , posicao inicial da sequencia , posicao final da sequencia
"""
def readNCBI(user,startSeq,stopSeq):
  Entrez.email = user
  handle = Entrez.efetch(db="nucleotide",id="NC_002942.5",rettype="gbwithparts",retmode="text",seq_start=startSeq,seq_stop=stopSeq)
  return handle

"""
 funcao para ir buscar mais informacoes de cada uma das locus tags
"""
def readProteinTag(user, genID):
  Entrez.email = user
  handle = Entrez.efetch(db="Genome",id=genID,rettype="gbwithparts",retmode="text")
  return handle

# escrever a sequencia para um ficheiro
def writeSequence(user,startSeq,stopSeq):
  handle = readNCBI(user,startSeq,stopSeq)
  record = SeqIO.read(handle,"genbank")
  SeqIO.write(record,"teste.gb","genbank")

"""
  escrever a informacao da locus tag para cada um dos ficheiros
"""
def writeProteinTag(user,genID,tag):
  handle = readProteinTag(user,genID)
  record = SeqIO.read(handle,"genbank")
  SeqIO.write(record,"locTags/"+tag+".gb","genbank")

"""
  Funcao trabalha com a sequencia.
  Le as features e guarda as informacoes que sao consideradas importantes - db_ref , function e locus_tag

  Falta talvez adicionar o ID da proteina e mais cenas
"""
def locateFeatures():
  records = SeqIO.read(open("example.gb"),"genbank")
  dic = {}
  flag = False
  for feature in records.features:
    #vale a pena considerar a source e o STS presente na nossa proteina
    if "gene" in feature.type and flag == False:
      qualifiers = feature.qualifiers
      flag = True
      if "locus_tag" in qualifiers:
        aux = Feature(qualifiers['locus_tag'][0])
        if "db_xref" in qualifiers.keys():
          aux.addDbRef(qualifiers['db_xref'][0])
        tag = qualifiers['locus_tag'][0]
        dic[tag] = aux
    elif flag == True:
      qualifiers = feature.qualifiers
      flag = False
      if "locus_tag" in qualifiers:
        aux = qualifiers['locus_tag'][0]
        if "function" in qualifiers.keys():
          dic[aux].addFunction(qualifiers['function'])
        if 'translation' in qualifiers.keys():
          dic[aux].addSeq(qualifiers['translation'][0])
        if 'note' in qualifiers.keys():
          dic[aux].addNote(qualifiers['note'][0])
        if 'protein_id' in qualifiers.keys():
          dic[aux].addProteinID(qualifiers['protein_id'][0])
        if 'product' in qualifiers:
          dic[aux].addProduct(qualifiers['product'][0])
    else:
      pass
  return dic

"""
  Funcao para obter valores adicionais de cada um dos genes -- ja nao me lembro onde xD
"""
def getInfoFeature(feature,tag):
  user  = "history.a61043@alunos.uminho.pt"
  if hasattr(feature,"id"):
    pass
    #writeProteinTag(user,feature.id,tag)
  else:
    #ESTA PARTE DO CODIGO NAO ESTA A SAKAR PARTES CORRECTAS - SAKA GENES DIFERENTES
    #COLOCAR UM PRINT E SABER QUAIS SAO E USAR A FUNCAO INDIVIDUAL PARA SAKAR OS GENES QUE FALTAM WriteInfoTag
    genID = re.split(r':',feature.db_xref)
    writeProteinTag(user,genID,tag)

"""
  Iterador de cada uma das proteinas na legionella
"""
def getProteinFeatures(dic):
  for feature in dic:
    getInfoFeature(dic[feature],feature)

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


"""
  funcao para realizar os blasts as sequencias dos locustag
"""
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

"""
  Passagem de informacao das features
"""
def stringInfo(feature):
  text = ""
  if hasattr(feature,"db_xref"):
    text += "db_xref: "+ feature.db_xref + "\n"
  if hasattr(feature,"function"):
    text += "Function: "+ ''.join(feature.function) + "\n"
  if hasattr(feature,"note"):
    text += "Note: " + feature.note + "\n"
  if hasattr(feature,"id"):
    text += "Protein_id: " + feature.id + "\n"
  if hasattr(feature,"product"):
    text += "Product: " + feature.product + "\n"
  text += " \n"
  return text

"""
  Comparar a ProteinTable da legionella com os resultados de todas as proteinas
  PODE
"""
def listProteinTable(dicto):
  with open("ProteinTable416_0.txt") as fp:
    flag = 0
    ofile = open("teste.txt","w")
    for line in fp:
      if flag > 0:
        text = line.split()
        aux = text[11:]
        genID = "GeneID:"+text[5]
        if genID in dicto:
          ofile.write(genID+"\n")
          ofile.write(' '.join(aux)+"\n")
          if hasattr(dicto[genID],'product'):
            func = ''.join(dicto[genID].product)
            ofile.write(func+"\n\n")
          else:
            ofile.write("\n")
      flag += 1

"""
  Funcao para sakar as informacoes dos genes protein id
"""
def DownLMissingGenID():
  dic = {"lpg2746":"19834312","lpg2729":"YP_096732"}
  user  = "history.a61043@alunos.uminho.pt"
  for arg in lista:
    writeInfoTag(user,arg,tag)

"""
  MAIN
"""
def main():
  folders = ["fastas","blasts","locTags"]
  createFolders(folders)

  #writeSequence("history.a61043@alunos.uminho.pt","2873801","3124550")

  dic = locateFeatures()
  getProteinFeatures(dic)
  #saveGenes(dic)
  #fastaToBlastFile("fastas/GeneID:19834107.fasta")
  #blastAllSeq("fastas")
  #getSwissProtInfo("19834166")
  #listProteinTable(dic)
  #readBlastProtein()
  #getUniProt()

main()
