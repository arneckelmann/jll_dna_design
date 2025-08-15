from Bio.SeqUtils import MeltingTemp as mt
from Bio import SeqIO
import xlsxwriter
from io import StringIO
import csv


import ipywidgets as widgets
from IPython.display import display



user_inputs = []
output = widgets.Output()



def ssm_ui():
    codon_positions_upload = widgets.FileUpload(
        accept='.csv',
        description='Position File',  # Accept all file types, or e.g. '.csv', '.txt'
        multiple=False  # Set to True to allow multiple files
    )


    genbank_file_upload = widgets.FileUpload(
        accept='.gbk',
        description='GenBank File',
        multiple=False
    )

    codon_text_box = widgets.Text(
        value='',
        placeholder='Mutagenesis Codon',
        description='Mutation',
        disabled=False
    )

    ORF_textbox = widgets.Text(
        value='',
        placeholder='CDS',
        description='Annotation',
        disabled=False
    )
    
    KLD_Gibson_option = widgets.RadioButtons(
        options=['KLD', 'Gibson'],
        description='Assembly Strategy:',
        disabled=False
    )
    

    primer_tm_text_box = widgets.FloatText(
        description='Tm:',
        disabled=False
    )
        
    submit_button = widgets.Button(
    value=False,
    description='Submit',
    disabled=False,
    button_style='', # 'success', 'info', 'warning', 'danger' or ''
    tooltip='Description',
    icon='check' # (FontAwesome names without the `fa-` prefix)
    )



    vbox = widgets.VBox([
                         genbank_file_upload, 
                         codon_positions_upload,
                         codon_text_box,
                         KLD_Gibson_option,
                         primer_tm_text_box, 
                         ORF_textbox,
                        submit_button],layout=widgets.Layout(
            width='50%', 
            border='solid 1px gray', 
            padding='10px', 
            align_items='stretch'
        ))
    
    
    display(vbox, output)
    
    user_inputs = [genbank_file_upload,
                   codon_positions_upload,
                  codon_text_box,
                   KLD_Gibson_option,
                   primer_tm_text_box,
                  ORF_textbox]
    
    def on_button_clicked(b):
        with output:
            execute_design(user_inputs)

    submit_button.on_click(on_button_clicked)

        
if __name__ == "__main__":
    UI_inputs()
    
def execute_design(user_inputs):
    
    genbank_info = user_inputs[0].value
    genbank_file_name = list(genbank_info.keys())[0]
    genbank_content = genbank_info[genbank_file_name]['content']
    genbank_bytes = genbank_content.decode('utf-8')
    genbank_str = StringIO(genbank_bytes)
    seq_record = list(SeqIO.parse(genbank_str, "genbank"))[0]
    features = seq_record.features
    
    
    label_coord = {}
    
    for feature in seq_record.features:
        if feature.type == 'CDS':
            label = feature.qualifiers['label'][0]
            dna_region = [list(feature.location)[0],list(feature.location)[-1]+1]
            label_coord[label] = dna_region
            
    codon_file_info = user_inputs[1].value
    codon_filename = list(codon_file_info.keys())[0]
    byte_position_str = codon_file_info[codon_filename]['content']
    positions = byte_position_str.decode('utf-8')
    aa_positions = positions.split("\r\n")
    aa_positions = [int(aa_position) for aa_position in aa_positions]
    nt_position_interval = [[(position*3)+50,(position*3)+53] for position in aa_positions] 
    
    
    
    codon_choice = user_inputs[2].value
    assembly_method = user_inputs[3].value
    target_tm = user_inputs[4].value
    
            
    input_annotation = user_inputs[-1].value
    start, end = label_coord[input_annotation]
    
    ORF_seq = str(seq_record.seq)[start-50:end+50]

    IDT_primer_table = design_mutagenesis_primers(ORF_seq,assembly_method, nt_position_interval, codon_choice, aa_positions, target_tm)
    generate_IDT_order(IDT_primer_table)
    
    
def rev_comp(sequence):
    sequence = sequence.upper()
    return sequence.replace("A","t").replace("T","a").replace("G","c").replace("C","g")[::-1].upper()

def get_primer_names(start,end,nt_position_interval,codon_choice, aa_positions):
    element_index = nt_position_interval.index([start,end])
    aa_position = aa_positions[element_index]
    primer_F_name = str(aa_position)+"_F"
    primer_R_name = str(aa_position)+"_R"
    return primer_F_name, primer_R_name


def design_mutagenesis_primers(ORF_seq,assembly_method, nt_position_interval, codon_choice, aa_positions, target_tm):
    
    f_primer_name_seq = []
    r_primer_name_seq = []
    
    for start, end in nt_position_interval:
        if assembly_method == "Gibson":
            interval_F1 = ORF_seq[start-24:start]

            interval_F2 = ORF_seq[end:end+17]
            base_count = 1

            while mt.Tm_NN(interval_F2) < target_tm - 2:
                interval_F2 = ORF_seq[end:end+17+base_count]
                base_count+=1

            full_F_primer_seq = interval_F1 + codon_choice + interval_F2
            full_R_primer_seq = rev_comp(interval_F1)

            base_trim_count = 1
            while mt.Tm_NN(full_R_primer_seq) > target_tm + 2:
                full_R_primer_seq = full_R_primer_seq[:base_trim_count*-1]
                base_trim_count += 1

            base_extension_count = 1
            while mt.Tm_NN(full_R_primer_seq) < target_tm - 2:
                interval_F1_seed = ORF_seq[start-24-base_extension_count:start]
                full_R_primer_seq = rev_comp(interval_F1_seed)
                base_extension_count += 1

            primer_F_name, primer_R_name = get_primer_names(start,end,nt_position_interval,codon_choice, aa_positions)
            
            

            f_primer_name_seq.append([primer_F_name,full_F_primer_seq])
            r_primer_name_seq.append([primer_R_name,full_R_primer_seq])


        elif assembly_method == "KLD":

            fwd_primer_seed = ORF_seq[end:end+17]
            rev_primer_seed = rev_comp(ORF_seq[start-17:start])

            fwd_extension_count = 1
            rev_extension_count = 1

            while mt.Tm_NN(fwd_primer_seed) < target_tm - 2:
                fwd_primer_seed = ORF_seq[end:end+17+fwd_extension_count]
                fwd_extension_count += 1

            while mt.Tm_NN(rev_primer_seed) < target_tm - 2:
                rev_primer_seed = rev_comp(ORF_seq[start-17-rev_extension_count:start])
                rev_extension_count += 1

            full_F_primer_seq = codon_choice + fwd_primer_seed
            full_R_primer_seq = rev_primer_seed

            primer_F_name, primer_R_name = get_primer_names(start,end,nt_position_interval,codon_choice)
            
            f_primer_name_seq.append([primer_F_name,full_F_primer_seq])
            r_primer_name_seq.append([primer_R_name,full_R_primer_seq])
            
    vbb_primer_name_seqs = design_vector_primers(ORF_seq,assembly_method)
    primer_name_seqs = f_primer_name_seq + r_primer_name_seq + vbb_primer_name_seqs
    
    primer_tables = [primer_name_seqs[i:i+96] for i in range(0,len(primer_name_seqs),96)]
    IDT_primer_table = {"Sheet"+str(i+1):[] for i in range(len(primer_tables))}
    wells_96 = [letter+str(i) for i in range(1,13) for letter in "ABCDEFGH"]
    header = ["Well Position","Name","Sequence"]
    
    element_count = 0
    
    for element in primer_tables:
        new_element = []
        row_count = 0
        for row in element:
            #row_index = element.index(row)
            current_well = wells_96[row_count]
            new_row = [current_well]+row
            new_element.append(new_row)
            row_count += 1
        element_count += 1
        IDT_primer_table["Sheet"+str(element_count)] = [header] + new_element
    
    return IDT_primer_table


def design_vector_primers(ORF_seq,assembly_method):
    
    if assembly_method == "Gibson":
        F1_interval = ORF_seq[26:50]
        F2_interval = ORF_seq[50:68]
        f_base_count = 1
        
        while mt.Tm_NN(F2_interval) > target_tm + 2:
            F2_interval = ORF_seq[50:68+f_base_count]
            f_base_count += 1
            
        VBB_F_primer = F1_interval+F2_interval
        
        R1_interval = rev_comp(ORF_seq[-68:-50])
        R2_interval = rev_comp(ORF_seq[-50:-24])
        
        r_base_count = 1
        while mt.Tm_NN(R1_interval) < target_tm - 2:
            R1_interval = rev_comp(ORF_seq[-68-r_base_count:-50])
            r_base_count += 1
        
        VBB_R_primer = R1_interval + R2_interval
        
        vbb_primer_name_seqs = [["VBB_F",VBB_F_primer],
                                ["VBB_R",VBB_R_primer]]
    else:
        vbb_primer_name_seqs = []
    
    
    
    return vbb_primer_name_seqs
    
def generate_IDT_order(IDT_primer_table):
    workbook = xlsxwriter.Workbook('output.xlsx')
    
    for sheet_name, sheet_data in IDT_primer_table.items():
        worksheet = workbook.add_worksheet(sheet_name)
        for row_idx, row in enumerate(sheet_data):
            for col_idx, cell in enumerate(row):
                worksheet.write(row_idx, col_idx, cell)
    
    
    workbook.close()