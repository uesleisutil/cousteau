# ********  DECLARACOES INICIAIS
import os
import Tkinter, tkFileDialog
import sys
import requests
import re
import shutil
from bs4 import BeautifulSoup

# By Vitor

# ABRE ARQUIVO DE ENTRADA
root    = Tkinter.Tk()
entrada = tkFileDialog.askopenfile(mode='r')
root.destroy()

#****************---------------correcao de bug--------------********************
if (entrada == None):
    sair = raw_input('\tArquivo de entrada nao selecionado. \n\t\tPressione enter para sair.\n')
    sys.exit()
#****************---------------fim da correcao--------------********************

pathname = os.path.dirname(entrada.name) #define o path de trabalho igual ao do arquivo de entrada
os.chdir(pathname)  #muda caminho de trabalho.

VALORES = []

# By Jean

while True:

    conteudo_linha = entrada.read().split("\n")
    VALORES.append(conteudo_linha)

    if (len(conteudo_linha) <= 1):
        break

print VALORES, "\n"


#### By Arthur

class Hidroweb(object):

    url_estacao = 'http://hidroweb.ana.gov.br/Estacao.asp?Codigo={0}&CriaArq=true&TipoArq={1}'
    url_arquivo = 'http://hidroweb.ana.gov.br/{0}'

    def __init__(self, estacoes):
        self.estacoes = estacoes

    def montar_url_estacao(self, estacao, tipo=1):
        return self.url_estacao.format(estacao, tipo)

    def montar_url_arquivo(self, caminho):
        return self.url_arquivo.format(caminho)

    def montar_nome_arquivo(self, estacao):
        return u'{0}.zip'.format(estacao)

    def salvar_arquivo_texto(self, estacao, link):
        r = requests.get(self.montar_url_arquivo(link), stream=True)
        if r.status_code == 200:
            with open(self.montar_nome_arquivo(estacao), 'wb') as f:
                r.raw.decode_content = True
                shutil.copyfileobj(r.raw, f)
            print '** %s ** (baixado)' % (estacao, )
        else:
            print '** %s ** (problema)' % (estacao, )

    def obter_link_arquivo(self, response):
        soup = BeautifulSoup(response.content)
        return soup.find('a', href=re.compile('^ARQ/'))['href']

    def executar(self):
        post_data = {'cboTipoReg': '10'}

        for est in self.estacoes:
            print '** %s **' % (est, )
            r = requests.post(self.montar_url_estacao(est), data=post_data)
            link = self.obter_link_arquivo(r)
            self.salvar_arquivo_texto(est, link)
            print '** %s ** (concluído)' % (est, )

if __name__ == '__main__':
    estacoes = VALORES[::1][0]
    hid = Hidroweb(estacoes)
    hid.executar()
