"""
  Parte BioPython
"""

from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import ExPASy
from Bio import SwissProt
from Bio.ExPASy import ScanProsite
from Bio.ExPASy import Prosite
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.ExPASy import Enzyme
from Bio.KEGG import REST

"""
  Imports genericos
"""
import os
import re
import time
import xml.etree.ElementTree as ET

"""
  class da feature
"""
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

  def addRefSeq(self,ref):
    self.ref = ref

  def addReview(self,rev):
    self.rev = rev

  def addDescpt(self,desc):
    self.desc = desc

  def geneName(self,genName):
    self.genName = genName

  def addKeyWords(self,keywords):
    self.keywords = keywords

  def addFeatures(self,features):
    self.features = features

  def addComments(self,comments):
    self.comments = comments

  def tamanhoSeq(self,tam):
    self.tam = tam

def locateFeatures():
  """
    Funcao trabalha com a sequencia.
    Le as features e guarda as informacoes que sao consideradas importantes - db_ref , function e locus_tag
  """
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


def getInfoFeature(feature,tag):
  """
    Funcao para obter valores adicionais de cada um dos genes -- ja nao me lembro onde xD
    esta funcao ja nao e necessaria
  """
  user  = "history.a61043@alunos.uminho.pt"
  if hasattr(feature,"id"):
    pass
    #writeProteinTag(user,feature.id,tag)
  else:
    #ESTA PARTE DO CODIGO NAO ESTA A SAKAR PARTES CORRECTAS - SAKA GENES DIFERENTES
    #COLOCAR UM PRINT E SABER QUAIS SAO E USAR A FUNCAO INDIVIDUAL PARA SAKAR OS GENES QUE FALTAM WriteInfoTag
    genID = re.split(r':',feature.db_xref)
    writeProteinTag(user,genID,tag)

######################################################################################################################

"""
  UNIPROT Functions
"""

def getUniProt(dic,tag):
  """
    obter informacao uniprot de uma unica prot - ver com as meninas
  """
  d = dic[tag]
  ref = d.ref
  handle = ExPASy.get_sprot_raw(ref)
  record = SwissProt.read(handle)

  d.addReview(record.data_class)
  d.addDescpt(record.description)
  d.geneName(record.gene_name)
  d.addKeyWords(record.keywords)
  d.addFeatures(record.features)
  d.addComments(record.comments)
  d.tamanhoSeq(record.sequence_length)

  """
  ['__doc__', '__init__', '__module__', 'accessions', 'annotation_update',
  'comments', 'created', 'cross_references', 'data_class', 'description',
  'entry_name', 'features', 'gene_name', 'host_organism', 'keywords',
  'molecule_type', 'organelle', 'organism', 'organism_classification',
  'references', 'seqinfo', 'sequence', 'sequence_length',
  'sequence_update', 'taxonomy_id']
  """

def getUniProtLegionella(dic):
  """
    Procurar tdas as proteinas da legionella na uniprot
  """
  for tag in dic:
    if hasattr(dic[tag],"ref"):
      if dic[tag].ref == "null":
        pass
    else:
      getUniProt(dic,tag)
      time.sleep(10)

######################################################################################################################

"""
  PDB Functions
"""


######################################################################################################################

"""
  BLAST FUNCTIONS
"""

"""
  funcao para realizar os blasts as sequencias dos locustag
"""
def blastAllSeq(folder):
  fastas = []
  fastas += [each for each in os.listdir(folder) if each.endswith('.fasta')]

  for fasta in fastas:
    fastaToBlastFile(fasta)
    time.sleep(300) #pausar a funcao por 3 minutos para nao estourar



def stringInfo(feature):
  """
    Passagem de informacao das features para
  """
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

######################################################################################################################

"""
  ITERATORS
"""

def writeGeneXML(user,dic):
  """
    Iterador para ir buscar o xml correspondente de cada gene da legionella
  """
  for feature in dic:
    genID = re.split(r':',dic[feature].db_xref)
    writeSequenceXML(user,genID[1],feature)

def getProteinFeatures(dic):
  """
  Iterador de cada uma das proteinas na legionella
  """
  for feature in dic:
    getInfoFeature(dic[feature],feature)

######################################################################################################################

"""
  LEITURAS
"""

def readNCBI(user,startSeq,stopSeq):
  """
  Ler a sequencia do NCBI da legionella

  parametros - utilizador , posicao inicial da sequencia , posicao final da sequencia
  """
  Entrez.email = user
  handle = Entrez.efetch(db="nucleotide",id="NC_002942.5",rettype="gbwithparts",retmode="text",seq_start=startSeq,seq_stop=stopSeq)
  return handle

def readNCBIGene(user,idG):
  """
    Le informacao gene
  """
  Entrez.email = user
  handle = Entrez.efetch(db="gene",id=idG,rettype="gb",retmode="xml")
  return handle


def readProteinTag(user, genID):
  """
    funcao para ir buscar mais informacoes de cada uma das locus tags
  """
  Entrez.email = user
  handle = Entrez.efetch(db="Genome",id=genID,rettype="gbwithparts",retmode="text")
  return handle

######################################################################################################################


"""
  Funcoes de mexer com o xml
"""

def readXMLGene(folder,dic):
  """
    xml do gene
  """
  xmls = []
  xmls += [each for each in os.listdir(folder) if each.endswith('.xml')]
  for xml in xmls:
    ofile = open(folder+"/"+xml,"r")
    tree = ET.parse(ofile)
    root = tree.getroot()
    place = root.findall("./Entrezgene/Entrezgene_comments/Gene-commentary/Gene-commentary_comment/Gene-commentary/Gene-commentary_products/Gene-commentary/Gene-commentary_products/Gene-commentary/Gene-commentary_comment/Gene-commentary/Gene-commentary_comment/Gene-commentary/Gene-commentary_source/Other-source/Other-source_src/Dbtag/Dbtag_tag/Object-id/Object-id_str")
    locustag = re.split(r'\.',xml)[0]
    if not place:
      dic[locustag].addRefSeq("null")
    else:
     for pl in place:
      dic[locustag].addRefSeq(pl.text)


def readXMLBlast(file):
  ### esta e so um teste
  result_handle = open("blasts/GeneID:19834107.xml")
  blast_record = NCBIXML.read(result_handle)
  E_VALUE_THRESH = 0.04 ##ThreshHold
  for alignment in blast_record.alignments:
     for hsp in alignment.hsps:
         if hsp.expect < E_VALUE_THRESH:
             print("****Alignment****")
             print("sequence:", alignment.title)
             print("length:", alignment.length)
             print("e value:", hsp.expect)
             print(hsp.query[0:75] + "...")
             print(hsp.match[0:75] + "...")
             print(hsp.sbjct[0:75] + "...")

def prepFichClustalW(dic,tag):
  genID = dic[tag].db_xref
  lista = []
  ofile = open("teste.fasta","w")
  result_handle = open("blasts/"+genID+".xml")
  blast_record = NCBIXML.read(result_handle)
  E_VALUE_THRESH = 0.04 ##ThreshHold
  for alignment in blast_record.alignments:
     for hsp in alignment.hsps:
         if hsp.expect < E_VALUE_THRESH:
            s = re.split(r'>',alignment.title)[0]
            filo = s[s.find("[")+1:s.find("]")]
            if filo not in lista:
              lista.append(filo)
              subject = filo.replace(' ','_')
              ofile.write(">"+subject+":\n")
              ofile.write(hsp.sbjct[0:])
              ofile.write("\n\n")
  ofile.close()

######################################################################################################################

"""
  ESCRITAS
"""
def writeSequence(user,startSeq,stopSeq):
  """
   escrever a sequencia total da legionella para um ficheiro
  """
  handle = readNCBI(user,startSeq,stopSeq)
  record = SeqIO.read(handle,"genbank")
  SeqIO.write(handle,"teste.gb","xml")

def writeSequenceXML(user,idG,tag):
  """
   escrever a sequencia total da legionella para um ficheiro
  """
  handle = readNCBIGene(user,idG)
  ofile = open("geneXML/"+tag+".xml","w")
  ofile.write(handle.read())
  ofile.close()


def writeProteinTag(user,genID,tag):
  """
  escrever a informacao da locus tag para cada um dos ficheiros
  """
  handle = readProteinTag(user,genID)
  record = SeqIO.read(handle,"genbank")
  SeqIO.write(record,"locTags/"+tag+".gb","genbank")

def proteinSequencesToFasta(dic):
  """
  Funcao que pega nas proteinas que tem uma sequencia associada e cria um ficheiro FASTA
  """
  for tag in dic:
    if hasattr(dic[tag], 'seq'):
      ofile = open("fastas/"+tag+".fasta", "w")
      ofile.write(">" + tag + "\n" +dic[tag].seq + "\n")
      #do not forget to close it
      ofile.close()
    else:
      pass

def fastaToBlastFile(fasta):
  """
    Pega no blast da sequencia e vai buscar o blast associado e guarda o
  """
  ofile = open("fastas/"+fasta,"r")
  print("Reading .... {}".format(ofile))
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



def makeGenText(dictionary):
  """
    Criar tabela
  """
  ofile = open("genes.txt","w")
  for gene in dictionary:
    txt = stringInfo(dictionary[gene])
    ofile.write(gene+"\n")
    ofile.write(txt)
  ofile.close()


######################################################################################################################

"""
  MISC
"""


def createFolders(directories):
  """
    Create folders
  """
  for directory in directories:
    if not os.path.exists(directory):
      os.makedirs(directory)



######################################################################################################################

"""
  Random testes
"""
def abc(dic):
  for i in dic:
      print("{} -> {}".format(i,dic[i].ref))

def readDominios(user):
  """
  Ler a sequencia do NCBI da legionella

  parametros - utilizador , posicao inicial da sequencia , posicao final da sequencia
  """
  Entrez.email = user
  handle = Entrez.esearch(db="cdd",id="YP_096627",rettype="xml",retmode="xml")
  print(handle.read())

def get_ECnumber():
    handle = open("enzyme.dat")
    records = Enzyme.parse(handle)
    ecnumbers = [record["DE"] for record in records]
    print(ecnumbers)

#####################################################################################################################

"""
  PROSITE
"""
def scanProsite():
  handle = ScanProsite.scan(seq="MKVKRNIIKELRVSRGWTQEHLAHLCDCSVRTIARAEKDGVNSMETITALAAVFGVEASQLMGMDKIEQLQQKTDEGERYLIRLHSGKMIVESFIESLAYKIDYEIPRNSKDADLISSVIQEIKDWGEIWGDLDVGEQIRASYRLNEIISEIEENGLVLFGLRTREKLNIFSQEIEGFLCNLFISYKDNERIVILGSDSTEE",lowscore=1)
  result = ScanProsite.read(handle)
  for i in range(len(result)):
    print(result[i])

def testProsite():
  handle = ExPASy.get_prosite_entry('PS50943')
  ofile = open("test.html","wb")
  ofile.write(handle.read())
  ofile.close()

"""
  CLUSTAL - testar pc da ana / daniela
"""
def clustalW():
  clustalomega_cline = ClustalOmegaCommandline(infile="teste.fasta", outfile="out.txt", verbose=True, auto=True, force=True)
  clustalomega_cline()

##########################################

def getEC():
  handle = open("enzyme.dat")
  records = Enzyme.parse(handle)
  ecnumbers = [record["ID"] for record in records]
  print(ecnumbers)

def Entrezparsing():
  ### FALTA O STRAND
  ofile = open("geneXML/lpg2549.xml")
  tree = ET.parse(ofile)
  root = tree.getroot()
  place = root.findall("./Entrezgene/Entrezgene_comments/Gene-commentary/Gene-commentary_comment/Gene-commentary/Gene-commentary_products/Gene-commentary/Gene-commentary_products/Gene-commentary/Gene-commentary_comment/Gene-commentary/Gene-commentary_comment/Gene-commentary/Gene-commentary_source/Other-source/Other-source_src/Dbtag/Dbtag_db")
  dbID = root.findall("./Entrezgene/Entrezgene_comments/Gene-commentary/Gene-commentary_comment/Gene-commentary/Gene-commentary_products/Gene-commentary/Gene-commentary_products/Gene-commentary/Gene-commentary_comment/Gene-commentary/Gene-commentary_comment/Gene-commentary/Gene-commentary_source/Other-source/Other-source_src/Dbtag/Dbtag_tag/Object-id/Object-id_id")
  lens = len(place)
  for i in range(lens):
    if place[i].text == "CDD":
      print("{} --> {}".format(place[i].text,dbID[i].text))


def cddTest(user,idC):
  Entrez.email = user
  handle = Entrez.esummary(db="cdd",id=idC)
  print(handle.read())
######################################################################################################################

"""
  ZONA DA ANA - KEEP OUT
"""

def infoFile(feature):
  """
  Passagem de informaçao da função das proteínas
  """
  txt=''
  if hasattr(feature,"function"):
        txt += "Function: "+ ''.join(feature.function) + "\n"
  return txt


def creatFile(dictionary, start, stop):
  """
  Cria um ficheiro com o gene e a sua função, se esta existir
  """
  file = open("gene_func.txt","w")
  #for gene in dictionary:
  for i in range(start,stop):
    try:
      txt = infoFile(dictionary["lpg"+str(i)])
      if txt!='':
          file.write("lpg"+str(i)+",")
          file.write(txt)
      else:
          file.write("lpg"+str(i)+"\n")
    except KeyError as e:
      print("nao tem o lpg{}".format(i))
  file.close()

def creatTable():
    """
    Cria uma tabela colorida no excel com as funções das proteínas
    """
    workbook = xls.Workbook('gene_func.xlsx')
    worksheet = workbook.add_worksheet()
    worksheet.set_column('A:A', 20)
    bold=workbook.add_format({'bold':True})
    i=0
    j=0
    for line in open("gene_func.txt","r"):
        if 'Function' in line:
            if 'Transport' in line:
                c='red'
            elif 'Metabolism' in line:
                c='yellow'
            elif 'Toxin production' in line:
                c='blue'
            elif 'Viral functions' in line:
                c='pink'
            elif 'cell division' in line:
                c='brown'
            elif 'secretion' in line:
                c='orange'
            elif 'DNA/RNA degradation' in line:
                c='green'
            elif 'Signal transduction' in line:
                c='purple'
            elif 'Translation' in line:
                c='navy'
            elif 'Transcription factors' in line:
                c='gray'
            elif 'Replication' in line:
                c='lime'
            elif 'Transcription' in line:
                c='cyan'
            elif 'Biodegradation' in line:
                c='magenta'
            elif 'Biosynthesis' in line:
                c='silver'
        else:
            c='white'
        color=workbook.add_format({'bg_color': c})
        worksheet.write(j, i, line[0:7], color)
        i+=1
        if i>15:
            i=0
            j+=1
    workbook.close()

#########################################

def swiss():
  handle = open("uniprot_sprot.dat")
  for record in SwissProt.parse(handle):
    print(record.description)
    time.sleep(5)

def get_genes():
    genes_met=[]
    genes_trans=[]
    genes_fact=[]
    for line in open("gene_func.txt","r"):
        if 'Metabolism' in line:
            genes_met.append(line[0:7])
        elif 'Transport' in line:
            genes_trans.append(line[0:7])
        elif 'Transcription factors' in line:
            genes_fact.append(line[0:7])
    return genes_met, genes_trans, genes_fact


def kegg():
  request = REST.kegg_get("lpn:lpg2547")
  print(request.read())

  #records = Enzyme.parse(request)
  #print(records.entry)

"""
  MAIN
"""

def main():

  #dic = locateFeatures()
  """
  readXMLGene("geneXML",dic)
  getUniProtLegionella(dic)
  makeGenText(dic)
  """

  #blastAllSeq("fastas")
  #prepFichClustalW(dic,"lpg2549")

  #writeGeneXML("history.a61043@alunos.uminho.pt",dic)
  #getProteinFeatures(dic)
  #saveGenes(dic)

  #readBlastProtein()
  #test("history.a61043@alunos.uminho.pt")

  #swiss()
  kegg()

  #scanProsite()
  #testProsite()

  #get_ECnumber()
  #clustalW()
  #Entrezparsing()

  #cddTest("history.a61043@alunos.uminho.pt",278590)


  """

  """
main()
