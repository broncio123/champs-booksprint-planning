MARKERS_latex = ['ch','sec','subsec','fig','tab']
out = []

def get_label_elements(label):
    try:
        marker, label_name = [x.strip() for x in label.split(':')]
        return marker, label_name
    except:
        return 'None'

def get_labels_tags_dic(nb, labels_pattern_latex, MARKERS_latex = MARKERS_latex):
    """Generate a dictionary of all LaTeX defined labels in Jupyter Notebook"""
    dic = {}
    for marker in MARKERS_latex:
        dic[marker] = {}
        COUNTER = 0
        #############################################
        # Search for LaTeX label/captions/reference commands within cells' content (source)
        #############################################
        for i in range(len(nb['cells'])):
            cell = nb['cells'][i]
            for j in range(len(cell['source'])):
                item = cell['source'][j]
                #############################################
                # Look for defined label in cell content
                #############################################
                labels_matched_raw = re.findall(labels_pattern_latex, item)
                if labels_matched_raw:
                    #############################################
                    # Add label names/numbers to dictionary according to marker
                    #############################################
                    label = labels_matched_raw[0]
                    if get_label_elements(label) != 'None':
                        if  get_label_elements(label)[0] == marker:
                            COUNTER += 1
                            label_name = get_label_elements(label)[1]
                            dic[marker][label_name] = str(COUNTER)
    return dic

##########################################################################################
def find_replace_citations_in_cell(cell, citations_pattern_latex):
    for j in range(len(cell['source'])):
        item = cell['source'][j]
        citations_matched_raw = re.findall(citations_pattern_latex, item)
        bibliography_matched_raw = re.findall(bibliography_pattern_latex,item)
        
        if citations_matched_raw:
            item_modified = replace_citation_syntax(citations_matched_raw, item)
            cell['source'][j] = item_modified
            
        elif bibliography_matched_raw:
            item_modified = replace_bibliography_syntax(bibliography_matched_raw, item)
            cell['source'][j] = item_modified
    return cell

def replace_citation_syntax(citations_matched_raw, item):
    item_modified = item
    for match in citations_matched_raw:
        bib_tags_raw = re.split(",", match)
        bib_tags_clean = [x.strip() for x in bib_tags_raw]
        #############################################
        # Replace LaTeX syntax with jekyll's in cell items
        #############################################
        jekyll_syntax_items = ["{% cite"]+bib_tags_clean+['--file '+bib_filename]+["%}"]
        separator = ' '
        chunk_new = separator.join(jekyll_syntax_items)
        chunk_original = '\cite{'+match+'}'
        item_modified = item_modified.replace(chunk_original, chunk_new)
    return item_modified

def replace_bibliography_syntax(bibliography_matched_raw, item):
    # Only one single bibliography line exists
    bib_tags_raw = re.split(",", bibliography_matched_raw[0])
    bib_tags_clean = [x.strip() for x in bib_tags_raw]
    separator = ' '
    jekyll_syntax_items = ["{% bibliography"]+[separator.join(['--file', x]) for x in bib_tags_clean]+["--cited %}"]
    chunk_new = separator.join(jekyll_syntax_items)
    chunk_original = r'\bibliography{'+bibliography_matched_raw[0]+'}'
    item_modified = item.replace(chunk_original, chunk_new)
    return item_modified

##########################################################################################
def find_replace_refs_in_cell(cell, refs_pattern_latex):
    for j in range(len(cell['source'])):
        item = cell['source'][j]
        refs_matched_raw = re.findall(refs_pattern_latex, item)
        if refs_matched_raw:
            item_modified = replace_ref_syntax(refs_matched_raw, item)
            cell['source'][j]  = item_modified
    return cell

def replace_ref_syntax(refs_matched_raw, item):
    item_modified = item
    for label in refs_matched_raw:
        if get_label_elements(label) != 'None':
            chunk_original = r'\ref{'+label+'}'
            marker, label_name = get_label_elements(label)
            label_number = labels_tags_dic[marker][label_name]
            chunk_new = '['+marker+':'+label_number+'](#'+label+')'
            item_modified = item.replace(chunk_original, chunk_new)
            item = item_modified
        else:
            pass
    return item

##########################################################################################

def find_labels_captions_in_cell(cell, labels_pattern_latex, captions_pattern_latex):
    """Find raw LaTeX labels/captions defined in cell. This assumes a single pair is only present."""
    dic = {}
    for j in range(len(cell['source'])):
        item = cell['source'][j]
        labels_matched_raw = re.findall(labels_pattern_latex, item)
        captions_matched_raw = re.findall(captions_pattern_latex, item)
        if labels_matched_raw:
            dic['label'] = labels_matched_raw[0]
        if captions_matched_raw:
            dic['caption'] = captions_matched_raw[0]
    return dic
 
def replace_labels_captions_in_cell(cell, labels_pattern_latex, captions_pattern_latex):
    dic = find_labels_captions_in_cell(cell, labels_pattern_latex, captions_pattern_latex)
    for j in range(len(cell['source'])):
        item = cell['source'][j]
        labels_matched_raw = re.findall(labels_pattern_latex, item)
        captions_matched_raw = re.findall(captions_pattern_latex, item)
        
        if labels_matched_raw:
            item_modified = replace_label_syntax(item, dic)
            cell['source'][j] = item_modified
            
        if captions_matched_raw:
            item_modified = replace_caption_syntax(item, dic)
            cell['source'][j] = item_modified
        else:
            item_modified = item
    return cell

def replace_label_syntax(item, dic):
    label_original = dic['label']
    try: 
        marker, label_name = get_label_elements(label_original)
        chunk_new = '\n<a id="'+label_original+'"></a>'
        item_modified = re.sub(labels_pattern_latex, chunk_new, item)
        return item_modified
    except:
        return item

def replace_caption_syntax(item, dic):
    caption = dic['caption']
    label_original = dic['label']
    marker, label_name = get_label_elements(label_original)
    label_number = labels_tags_dic[marker][label_name]
    #############################################
    html_syntax_items = (
        '<figcaption style="text-align:center;font-size:14px">',
        '<b>',
        marker+':'+label_number+' ',
        '</b>',
        '<em>'+' ',
        caption,
        '</em>',
        '</figcaption>'
    )
    chunk_new = ''.join(html_syntax_items)
    chunk_original = r'\caption{'+caption+'}'
    item_modified = item.replace(chunk_original, chunk_new)
    return item_modified
##########################################################################################

def latex_to_html_in_jupyter(nb):
    for i in range(len(nb['cells'])):
        cell = nb['cells'][i]
        cell_modified = find_replace_citations_in_cell(cell, citations_pattern_latex)
        cell_modified = find_replace_refs_in_cell(cell_modified, refs_pattern_latex)
        cell_modified = replace_labels_captions_in_cell(cell, labels_pattern_latex, captions_pattern_latex)
        nb['cells'][i] = cell_modified
    return nb

##########################################################################################
if __name__ == "__main__":
    import re
    import sys
    import json
    #############################################
    # Settings
    #############################################
    citations_pattern_latex = re.compile(r'\\cite\{(.+?)\}')
    bibliography_pattern_latex = re.compile(r'\\bibliography\{(.+?)\}')
    refs_pattern_latex = re.compile(r'\\ref\{(.+?)\}')
    labels_pattern_latex = re.compile(r'\\label\{(.+?)\}$')
    captions_pattern_latex = re.compile(r'\\caption\{(.+?)\}$')
    #############################################
    # Script input arguments (self-explanatory)
    #############################################
    nb_infile = sys.argv[1]
    nb_outfile = sys.argv[2]
    bib_filename = sys.argv[3]
    #############################################
    # Load Jupyter Notebook as Dictionary
    #############################################
    try:
        with open(nb_infile,'r') as fp:
            nb = json.load(fp)
        fp.close()
        #############################################
        # Turn LaTeX label/ref syntax in Jupyter Noteboks into MD syntax
        #############################################
        labels_tags_dic = get_labels_tags_dic(nb, labels_pattern_latex)
        nb_modified =  latex_to_html_in_jupyter(nb)
        
        with open(nb_outfile,'w') as fp:
            json.dump(nb_modified, fp)
        fp.close()
        print("New Notebook successfully generated")
    except:
        print("Couldn't find Jupyter Notebook. Check your input path")
