#!/usr/bin/python

from Bio import SeqIO
import numpy as np
import pandas as pd
import os.path
from os import path

basedir='/scratch/brown/ggavelis/all_dinos/'
full_ips_dir='/scratch/brown/ggavelis/all_dinos/AllInterProScan/full_tsv/'. #Directory where input '.tsv' is located
summary_ips_dir='/scratch/brown/ggavelis/all_dinos/AllInterProScan/summary_tsv/' #Directory for output

def set_to_commalist(example):  #This is gonna clean up the formatting of sets before I put them in my output tsv as comma-sep lists
    return str(example).replace('{\'\', ','').replace('{','').replace('}','').replace('\'','').replace(', ,',',')

verbiage=[' profile', ' signature', ' family', ' Family', ' domain', ' superfamily', '.', '-like domain superfamily', ',', '-like', 'Protein of unknown function ', 'Domain of unknown function ', ' N-terminal', ' C-terminal', 'Putative ', '-containing', 'EF-hand calcium-binding', 'EF-Hand 1 calcium-binding site', 'EF-hand pair protein CML', ' type', '-type', 'EF-hand pair'] 

summary_columns=[]
ips_columns=['accession', 'MD5', 'length', 'refType', 'FeatureAccession', 'FeatureDescription', 'start_loc', 'stop_loc', 'score', 'status', 'date', 'iprlookup_accesion', 'iprlookup_description', 'goterms', 'pathways']
#summary_columns=['query', 'IPS_annotations', 'transmembrane', 'GO_terms', 'pathways', 'IPS_annot_types', 'IPS_annot_accessions']
summary_columns='query'+'\t'+'IPS_annotations'+'\t'+'transmembrane'+'\t'+'GO_terms'+'\t'+'pathways'+'\t'+'IPS_annot_types'+'\t'+'IPS_annot_accessions'
basedir='/scratch/brown/ggavelis/all_dinos/'

libraries = [ \
'Dinophysis_norvegica_DN2_GG~120300~2864~', \
'Dinophysis_norvegica_DN1_GG~120300~2864~', \
'Dinophysis_norvegica_DN4_GG~120300~2864~', \
'Dinophysis_norvegica_DN5_GG~120300~2864~', \
'Alexandrium_catenella_MMETSP~2925~2864~MMETSP0790', \
'Alexandrium_monilatum_MMETSP4~311494~2864~MMETSP0093,MMETSP0095,MMETSP0096,MMETSP0097', \
'Alexandrium_tamarense_MMETSP4~2926~2864~MMETSP0378,MMETSP0380,MMETSP0382,MMETSP0384', \
'Amphidinium_carterae_MMETSP4~2961~2864~MMETSP0259,MMETSP0398,MMETSP0399,MMETSP0267', \
'Amphidinium_massartii_MMETSP~160604~2864~MMETSP0689', \
'Amyloodinium_ocellatum0_SRA~79898~2864~', \
'Ansanella_granifera0_SRA~1872691~2864~', \
'Apocalathium_aciculiferum_MMETSP2~268820~2864~MMETSP0370,MMETSP0371', \
'Azadinium_spinosum_MMETSP3~632150~2864~MMETSP1036,MMETSP1037,MMETSP1038', \
'Brandtodinium_nutricula0_SRA~1333877~2864~', \
'Chromera_velia_MMETSP~505693~33630~MMETSP0290', \
'Crypthecodinium_cohnii_MMETSP4~2866~2864~MMETSP0323,MMETSP0324,MMETSP0325,MMETSP0326', \
'Cryptoperidiniopsis_sp0_GG~79893~2864~', \
'Dinophysis_acuminata_MMETSP~47934~2864~MMETSP0797', \
'Dinophysis_fortii0_SRA~150623~2864~', \
'Durinskia_baltica_MMETSP2~400756~2864~MMETSP0116,MMETSP0117', \
'Gambierdiscus_australes_MMETSP~439317~2864~MMETSP0766', \
'Gonyaulax_spinifera_MMETSP~66791~2864~MMETSP1439', \
'Green_Dinoflagellate_M0_SRA~2608282~2864~', \
'Green_Dinoflagellate_T0_SRA~2608281~2864~', \
'Gymnodinium_catenatum_MMETSP~39447~2864~MMETSP0784', \
'Gymnoxanthella_radiolariae0_SRA~1798043~2864~', \
'Gyrodiniellum_shiwhaense0_SRA~926437~2864~', \
'Heterocapsa_arctica_MMETSP~192219~2864~MMETSP1441', \
'Heterocapsa_rotundata_MMETSP~89963~2864~MMETSP0503', \
'Karenia_brevis_MMETSP2_U~156230~2864~MMETSP0201,MMETSP0202', \
'Karlodinium_veneficum_MMETSP3~407301~2864~MMETSP1015,MMETSP1016,MMETSP1017', \
'Kryptoperidinium_foliaceum_MMETSP2_CCMP~160619~2864~MMETSP0120,MMETSP0121', \
'Lepidodinium_chlorophorum0_SRA~107758~2864~', \
'Lingulodinium_polyedra_MMETSP4~160621~2864~MMETSP1032,MMETSP1033,MMETSP1034,MMETSP1035', \
'Margalefidinium_polykrikoides0_SRA~77300~2864~', \
'Noctiluca_scintillans_MMETSP~2966~2864~MMETSP0258', \
'Nusuttodinium_aeruginosum0_SRA~1488663~2864~', \
'Oxyrrhis_marina0_MM~2969~33630~', \
'Paragymnodinium_shiwhaense0_SRA~409284~2864~', \
'Pelagodinium_beii0_SRA~43686~2864~', \
'Peridinium_bipes0_SRA~2868~2864~', \
'Pfiesteria_piscicida0_SRA~71001~2864~', \
'Polykrikos_lebouriae0_GG~370573~2864~', \
'Proterythropsis_metatranscriptome0_GG~651421~2864~', \
'Prorocentrum_minimum_MMETSP3~39449~2864~MMETSP0268,MMETSP0269,MMETSP0270', \
'Pyrocystis_lunula2_SRA~2972~2864~', \
'Ross_Sea_Dinoflagellate0_SRA~2069452~2864~', \
'Scrippsiella_hangoei_MMETSP5~268821~2864~MMETSP0359,MMETSP0360,MMETSP0361,MMETSP0368,MMETSP0369', \
'Scrippsiella_trochoidea_MMETSP3~71861~2864~MMETSP0271,MMETSP0272,MMETSP0286', \
'Togula_jolla_MMETSP~285029~2864~MMETSP0224', \
'Yihiella_yeosuensis0_SRA~1744892~2864~', \
'Breviolum_minutum_vB~2499525~2864~', \
'Amoebophrya_sp_EP~88551~33630~', \
'Cladocopium_goreaui_vB~2562237~2864~', \
'Cladocopium_sp_C92_vB~2486705~2864~', \
'Fugacium_kawagutii_vB~2697096~2864~', \
'Perkinsus_chesapeaki_v1~330153~33630~', \
'Perkinsus_marinus_v1~31276~33630~', \
'Perkinsus_olseni_v1~32597~33630~', \
'Polarella_glacialis_CCMP1383_vB~89957~2864~', \
'Polarella_glacialis_CCMP2088_vB~89957~2864~', \
'Symbiodinium_microadriaticum_vB~2951~2864~', \
'Symbiodinium_tridacnidorum_vb~1602974~2864~', \
'Vitrella_brassicaformis_v1~1169539~2864~', \
]

### UNUSED TRANSCRIPTOMES
# 'Karenia_brevis_MMETSP2_Wilson~156230~2864', \ #Redundant w/ other karenia sample, but lower BUSCO coverage
# 'Alexandrium_fundyense_MMETSP3~2932~2864', \ #Systematically contaminated by diverse organisms
# 'Alexandrium_minutum_MMETSP~39455~2864', \  #Insufficient BUSCO
# 'Kryptoperidinium_foliaceum_MMETSP2_CCAP~160619~2864', \ #Redundant w/ conspecific, but this one has less BUSCO coverage

for i in libraries:
    
    
    count,ips_lines,annot_count,previous_accession = 0,0,0,''
    library,recipient,ancestral,mmetsp = i.split('~')
    print('\n***** ' + library + ' ******\n')
    
    out_file = summary_ips_dir + library + '.ips_summary.tsv'
    if not path.exists(out_file):
        out_handle = open(out_file, "w")
        out_handle.write(summary_columns+'\n')
    
    
        ### Find InterProScan
        ips_file = full_ips_dir + library + '_ips.tsv'
        if not path.exists(ips_file):
            print('no InterProScan tsv file for ' + library)
        
        ips_df_nan = pd.read_csv(ips_file, delimiter='\t', names=ips_columns)
        ips_df_raw = ips_df_nan.replace(np.nan, '', regex=True)                  # replace 'NaN' values with empty strings

        exclude_annot = ['MobiDBLite', 'TMHMM']    
        annot_count,ips_lines,TMcount = 0,0,0
        previous_accession=''
    
        ######## Annotation Info that I want to keep
        feature_description_set = set()
        ipr_lookup_description_set = set()
        goterms_set = set()
        pathways_set = set()
        reference_type_set = set()
        feature_accession_set = set()
        # iprlookup_accession_set = set()  #not useful
        transmembrane_set = set()

        ####### Go through raw IPS file
        for index, row in ips_df_raw.iterrows():
            ips_lines += 1
            query = row['accession']

            ###### Are we on the next protein yet?
            if row['accession'] != previous_accession:
                previous_accession = (row['accession'])
                annot_count += 1
                if ' ' in feature_description_set:
                    #print(' removed empty feature from set')
                    feature_description_set.remove(' ')
            
                prev_row = query+'\t'+set_to_commalist(feature_description_set)+'\t'+set_to_commalist(transmembrane_set)+'\t'+set_to_commalist(goterms_set)+'\t'+set_to_commalist(pathways_set)+'\t'+set_to_commalist(reference_type_set)+'\t'+set_to_commalist(feature_accession_set)#+'\t'+set_to_commalist(iprlookup_accession_set)
                prev_row = prev_row.replace('set()','') #remove empty sets
                prev_row = prev_row.replace(',\t','\t') #remove empty list elements
                prev_row = prev_row.replace(', \t','\t') #remove blank list elements
        
                ### If so, refresh variables
                feature_description_set = set()
                ipr_lookup_description_set = set()
                goterms_set = set()
                pathways_set = set()
                reference_type_set = set()
                feature_accession_set = set()
                iprlookup_accession_set = set()
                transmembrane_set = set()
        
                if ips_lines != 1:
                    out_handle.write(prev_row+'\n')
                    #print(prev_row+'\n')
    
            ###### By storing annotation info in sets (rather than lists), we remove redundancies (e.g. Pfam and PANTHER might both call the feature 'receptor kinase')
            if row['goterms'] != '':
                go_termlist = row['goterms'].split('|')   ### reformat '|' delimited Goterms into proper ',' delimited list
                for go in go_termlist:
                    goterms_set.add(go)  ### Now collapse the listed terms into a set
    
            if row['pathways'] != '':
                pathways_list = row['pathways'].split('|') ## reformat '|' delimited pathways into proper ',' delimited list
                for pathway in pathways_list:
                    if len(pathways_set) < 4:
                        pathways_set.add(pathway)   #Limit 3 pathways maximum to avoid over-annotating this (e.g. proteasomes have dozens of reactomes. Don't care)
                feature_accession_set.add(row['FeatureAccession'])
                reference_type_set.add(row['refType'])
        
            if row['FeatureAccession'] == 'TMhelix':
                transmembrane_set.add('Y')
                reference_type_set.add('TMHMM')
                TMcount+=1
        
            if row['refType'] not in exclude_annot:                  #### remove redundant verbiage from feature descriptions
                feature_description = row['FeatureDescription']
                feature_description = feature_description.replace('Ribosomal protein ','Ribosomal ').replace(' conserved site',' site').replace(' enzymes',' enzyme').replace(' region',' site')
                for word in verbiage:
                    feature_description = feature_description.replace(word,'')
                feature_description_set.add(feature_description)
            
                feature_description = row['iprlookup_description']
                feature_description = feature_description.replace('Ribosomal protein ','Ribosomal ').replace(' conserved site',' site').replace(' enzymes',' enzyme').replace(' region',' site')
                for word in verbiage:
                    feature_description = feature_description.replace(word,'')
                feature_description_set.add(feature_description)            
            
                if feature_description != '':   #### If this annotation method gave us a feature description, save accession and metadata
                    reference_type_set.add(row['refType'])
                    feature_accession_set.add(row['FeatureAccession'])
                    iprlookup_accession_set.add(row['iprlookup_accesion'])
            
                feature_description_set.add(feature_description)

        
        
        print(str(annot_count) + ' proteins annotated')
        print(str(TMcount) + ' had transmembrane helices')
        print('read across ' + str(ips_lines) + ' lines')
        print('written to ' + out_file)
        out_handle.close()
        
    else:
        print(library + ' already has an outfile')
        
print('Process complete')
