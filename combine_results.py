
import matplotlib
matplotlib.use('Agg')
import os, json, setlog
import matplotlib.pyplot as plt
import pickle
import ast
global combined_log
combined_log = setlog.init(toconsole=True, logname="combinedlog")

def combine_core_results(result_dict, outdir):
    combined_log.info("Core result combination start")
    if os.path.exists(outdir) != True:
        os.mkdir(outdir)
    if os.path.exists(os.path.join(outdir, "tables")) != True:
        os.mkdir(os.path.join(outdir, "tables"))
    json_result_file = os.path.join(os.path.join(outdir, "tables"), "combined_core_table.json")
    combined_table_file = open(os.path.join(os.path.join(outdir, "tables"),  "combined_core_table.tsv"), "w")
    combined_dict = {}
    attr_dict = {}
    identifier_line = ""
    core_gene_to_organism = {}
    total_org_count = float(len(result_dict.keys()))
    for res_dir in result_dict:
        input = result_dict[res_dir]
        tables_directory_path = os.path.join(res_dir, "tables")
        coretable_filepath = os.path.join(tables_directory_path, "coretable.tsv")
        core_file = open(coretable_filepath, "r")
        for line in core_file:
            if line.startswith("#"):
                identifier_line = line.strip("\n")
                continue
            else:
                split_line = line.split("\t")
                core_gene_with_desc = "\t".join(split_line[0:3])
                hit_dict = {"No" : "0", "Yes" : "1"}
                if core_gene_with_desc not in core_gene_to_organism:
                    core_gene_to_organism[core_gene_with_desc] = [input]
                else:
                    if input not in core_gene_to_organism[core_gene_with_desc]:
                        core_gene_to_organism[core_gene_with_desc].append(input)
                if core_gene_with_desc not in combined_dict:
                    attrs = []
                    hit_type = 0
                    for hit in split_line[3:7]:
                        if hit.lower() == "n/a":
                            hit = "No"
                        attrs.append(hit_dict[hit])
                        if hit_type == 0:
                            if hit == "Yes":
                                attr_dict[core_gene_with_desc + "|dup"] = [input]
                            else:
                                attr_dict[core_gene_with_desc + "|dup"] = []
                        if hit_type == 1:
                            if hit == "Yes":
                                attr_dict[core_gene_with_desc + "|bgc"] = [input]
                            else:
                                attr_dict[core_gene_with_desc + "|bgc"] = []
                        if hit_type == 2:
                            if hit == "Yes":
                                attr_dict[core_gene_with_desc + "|phl"] = [input]
                            else:
                                attr_dict[core_gene_with_desc + "|phl"] = []
                        if hit_type == 3:
                            if hit == "Yes":
                                attr_dict[core_gene_with_desc + "|knw"] = [input]
                            else:
                                attr_dict[core_gene_with_desc + "|knw"] = []
                        hit_type += 1
                    dup = str(attr_dict[core_gene_with_desc + "|dup"])
                    bgc = str(attr_dict[core_gene_with_desc + "|bgc"])
                    phl = str(attr_dict[core_gene_with_desc + "|phl"])
                    knw = str(attr_dict[core_gene_with_desc + "|knw"])
                    combined_dict[core_gene_with_desc] = core_gene_with_desc + "\t" + "\t".join(attrs) +"\t" + split_line[7].strip("\n") + "\t" +str(core_gene_to_organism[core_gene_with_desc]) +\
                                                         "\t" + dup + "\t"+ bgc + "\t"+ phl + "\t"+ knw + "\n"
                else:
                    attrs = combined_dict[core_gene_with_desc].split("\t")[3:7]
                    attr_number = 0
                    hit_type = 0
                    for hit in split_line[3:7]:
                        if hit.lower() == "n/a":
                            hit = "No"
                        attrs[attr_number] = str(int(attrs[attr_number]) + int(hit_dict[hit]))
                        attr_number += 1
                        if hit_type == 0:
                            if hit == "Yes":
                                attr_dict[core_gene_with_desc + "|dup"].append(input)
                        if hit_type == 1:
                            if hit == "Yes":
                                attr_dict[core_gene_with_desc + "|bgc"].append(input)
                        if hit_type == 2:
                            if hit == "Yes":
                                attr_dict[core_gene_with_desc + "|phl"].append(input)
                        if hit_type == 3:
                            if hit == "Yes":
                                attr_dict[core_gene_with_desc + "|knw"].append(input)
                        hit_type += 1
                    dup = str(attr_dict[core_gene_with_desc + "|dup"])
                    bgc = str(attr_dict[core_gene_with_desc + "|bgc"])
                    phl = str(attr_dict[core_gene_with_desc + "|phl"])
                    knw = str(attr_dict[core_gene_with_desc + "|knw"])
                    listed_hits = combined_dict[core_gene_with_desc].split("\t")[7]
                    combined_dict[core_gene_with_desc] = core_gene_with_desc + "\t" + "\t".join(attrs) + "\t" +listed_hits[0:-1] +"; "  + split_line[7][1:].strip("\n") + "\t" + str(core_gene_to_organism[core_gene_with_desc]) +\
                                                         "\t" + dup + "\t" + bgc + "\t" + phl + "\t" + knw  + "\n"
    identifier_line = identifier_line + "\t" + "Core_orgs" + "\t" + "Dup_orgs" + "\t" + "BGC_orgs" + "\t" + "Phyl_orgs" + "\t" + "Known_orgs" + "\n"
    combined_table_file.write(identifier_line)
    json_dict = {}
    for i in sorted(combined_dict):
        splitlist = combined_dict[i].split("\t")
        splitlist[3] = str(round(float(combined_dict[i].split("\t")[3]) / total_org_count, 2))
        splitlist[4] = str(round(float(combined_dict[i].split("\t")[4]) / total_org_count, 2))
        splitlist[5] = str(round(float(combined_dict[i].split("\t")[5]) / total_org_count, 2))
        splitlist[6] = str(round(float(combined_dict[i].split("\t")[6]) / total_org_count, 2))
        combined_dict[i] = "\t".join(splitlist)
        combined_table_file.write(combined_dict[i])
        json_dict[combined_dict[i].split("\t")[0]] = combined_dict[i].strip("\n").split("\t")[1:]
    combined_table_file.close()
    with open(json_result_file, 'w') as outfile:
        json.dump(json_dict, outfile, indent = 2)
    combined_log.info("Core result combination end")


def combine_known_results(result_dict, outdir):
    combined_log.info("Known result combination start")
    if os.path.exists(outdir) != True:
        os.mkdir(outdir)
    if os.path.exists(os.path.join(outdir, "tables")) != True:
        os.mkdir(os.path.join(outdir, "tables"))
    json_result_file = os.path.join(os.path.join(outdir, "tables"),  "combined_known_table.json")
    combined_table_file = open(os.path.join(os.path.join(outdir, "tables"),  "combined_known_table.tsv"), "w")
    known_hit_to_organism = {}
    known_hit_to_res_dir = {}
    for res_dir in result_dict:
        input = result_dict[res_dir]
        tables_directory_path = os.path.join(res_dir, "tables")
        knowntable_filepath = os.path.join(tables_directory_path, "knownhits.tsv")
        if os.path.exists(knowntable_filepath) != True:
            combined_log.debug("Knownhits table not found, exitting")
            return
        known_file = open(knowntable_filepath, "r")
        res_list = []
        for line in known_file:
            if line.startswith("#"):
                continue
            split_line = line.split("\t")
            res_model = split_line[0] + "\t" + split_line[1]
            if res_model in res_list:
                continue
            else:
                if res_model not in known_hit_to_organism:
                    known_hit_to_organism[res_model] = [input]
                else:
                    known_hit_to_organism[res_model].append(input)
                if res_model not in known_hit_to_res_dir:
                    known_hit_to_res_dir[res_model] = [res_dir]
                else:
                    known_hit_to_res_dir[res_model].append(res_dir)
                res_list.append(res_model)
    identifier_line = "#Model\tDescription\tOrg_count\tOrganisms\tPaths\n"
    combined_table_file.write(identifier_line)
    json_dict = {}
    for i in sorted(known_hit_to_organism):
        json_dict[i.split("\t")[0]] = [i.split("\t")[1]]
        json_dict[i.split("\t")[0]].append(str(len(known_hit_to_organism[i])))
        json_dict[i.split("\t")[0]].append(str(known_hit_to_organism[i]))
        json_dict[i.split("\t")[0]].append(str(known_hit_to_res_dir[i]))
        combined_table_file.write(i + "\t" + str(len(known_hit_to_organism[i])) + "\t")
        org = str(known_hit_to_organism[i])
        path = str(known_hit_to_res_dir[i])
        combined_table_file.write(org + "\t" + path + "\n")
    combined_table_file.close()
    with open(json_result_file, 'w') as outfile:
        json.dump(json_dict, outfile, indent = 2)
    combined_log.info("Known result combination end")

def combine_dup_results(result_dict, outdir):
    combined_log.info("Dup result combination start")
    if os.path.exists(outdir) != True:
        os.mkdir(outdir)
    if os.path.exists(os.path.join(outdir, "tables")) != True:
        os.mkdir(os.path.join(outdir, "tables"))
    json_result_file = os.path.join(os.path.join(outdir, "tables"), "combined_dup_table.json")
    combined_table_file = open(os.path.join(os.path.join(outdir, "tables"),  "combined_dup_table.tsv"), "w")
    core_hit_to_organism = {}
    core_hit_to_res_dir = {}
    for res_dir in result_dict:
        input = result_dict[res_dir]
        tables_directory_path = os.path.join(res_dir, "tables")
        duptable_filepath = os.path.join(tables_directory_path, "duptable.tsv")
        dup_file = open(duptable_filepath, "r")
        for line in dup_file:
            if line.startswith("#"):
                continue
            split_line = line.split("\t")
            core_gene = split_line[0] + "\t" + split_line[-1].strip("\n")
            if core_gene not in core_hit_to_organism:
                core_hit_to_organism[core_gene] = [input]
            else:
                core_hit_to_organism[core_gene].append(input)
            if core_gene not in core_hit_to_res_dir:
                core_hit_to_res_dir[core_gene] = [res_dir]
            else:
                core_hit_to_res_dir[core_gene].append(res_dir)
    identifier_line = "#Core_gene\tDescription\tOrg_count\tOrganisms\tPaths\n"
    combined_table_file.write(identifier_line)
    json_dict = {}
    for i in sorted(core_hit_to_organism):
        json_dict[i.split("\t")[0]] = [i.split("\t")[1]]
        json_dict[i.split("\t")[0]].append(str(len(core_hit_to_organism[i])))
        json_dict[i.split("\t")[0]].append(str(core_hit_to_organism[i]))
        json_dict[i.split("\t")[0]].append(str(core_hit_to_res_dir[i]))
        combined_table_file.write(i + "\t" + str(len(core_hit_to_organism[i])) + "\t")
        combined_table_file.write(str(core_hit_to_organism[i]) + "\t" + str(core_hit_to_res_dir[i]) + "\n")
    combined_table_file.close()
    with open(json_result_file, 'w') as outfile:
        json.dump(json_dict, outfile, indent = 2)
    combined_log.info("Dup result combination end")


def generate_summary(result_dict, outdir):
    if os.path.exists(outdir) != True:
        os.mkdir(outdir)
    if os.path.exists(os.path.join(outdir, "tables")) != True:
        os.mkdir(os.path.join(outdir, "tables"))
    json_result_file = os.path.join(os.path.join(outdir, "tables"), "summary_table.json")
    combined_table_file = open(os.path.join(os.path.join(outdir, "tables"),  "summary_table.tsv"), "w")
    organism_stats_dict = {}
    for res_dir in result_dict:
        total_genes = ""
        core_genes = ""
        total_bgc_hits = ""
        known_resistance = "N/A"
        gene_dup = "N/A"
        bgc_prox = "N/A"
        phyl = "N/A"
        two = ""
        three = ""
        input = result_dict[res_dir]
        log_file = open(os.path.join(res_dir, "arts-query.log"), "r")
        try:
            for line in log_file:
                if "CDS features: " in line:
                    total_genes = line.split("CDS features: ")[1].split(";")[0]
                    total_bgc_hits = line.split("Clusters: ")[1].strip()
                if "Known Resistance Hits: " in line:
                    known_resistance = line.split("Known Resistance Hits: ")[1].strip()
                if "DEBUG - extractdbgenes" in line and core_genes == "":
                    core_genes = line.split("):")[0].split()[-1]
                if "Found" in line and "duplicate genes" in line:
                    gene_dup = line.split("Found ")[1].split()[0]
                if "Proximity hits found: " in line:
                    bgc_prox = line.split("Proximity hits found: ")[1].strip()
                if "Phylogeny hits found: " in line:
                    phyl = line.split("Phylogeny hits found: ")[1].strip()
                if "Hits with two or more criteria: " in line:
                    two = line.split("Hits with two or more criteria: ")[1].split()[0]
                if "Hits with three or more criteria: " in line:
                    three = line.split("Hits with three or more criteria: ")[1].split()[0]
                organism_stats_dict[input] = [os.path.basename(res_dir),total_genes, total_bgc_hits, known_resistance, core_genes, gene_dup, bgc_prox, phyl, two, three]
            combined_log.info(os.path.join(res_dir, "arts-query.log") + " parsed succesfully")
        except Exception as e:
            combined_log.error("Something wrong with arts-query.log")
    identifier_line = "#Organism\tJobid\tTotal Genes\tTotal BGC Hits\tKnown Resistance\tCore Genes\tGene Duplication\tBGC Proximity\tPhylogeny / HGT\t2+\t3+\n"
    combined_table_file.write(identifier_line)
    for org in organism_stats_dict:
        combined_table_file.write(os.path.basename(org))
        for i in organism_stats_dict[org]:
            combined_table_file.write("\t")
            combined_table_file.write(i)
        combined_table_file.write("\n")
    combined_table_file.close()
    with open(json_result_file, 'w') as outfile:
        json.dump(organism_stats_dict, outfile, indent=2)

def generate_plots(maindir):
    bgc_names = []
    bgc_counts = []
    dup_names = []
    dup_counts = []
    core_names = []
    core_counts = []
    known_counts = []
    known_names = []
    bgc_statistics = {}
    dup_statistics = {}
    core_statistics = {}
    known_statistics = {}
    tables_dir = os.path.join(maindir, "tables")
    bgc_path = os.path.join(tables_dir, "combined_bgc_table.tsv")
    bgc_file = open(bgc_path, "r")
    core_path = os.path.join(tables_dir, "combined_core_table.tsv")
    core_file = open(core_path, "r")
    known_path = os.path.join(tables_dir, "combined_known_table.tsv")
    known_file = open(known_path, "r")
    dup_path = os.path.join(tables_dir, "combined_dup_table.tsv")
    dup_file = open(dup_path, "r")
    if os.path.exists(known_path):
        for line in known_file:
            if line.startswith("#"):
                continue
            line = line.split("\t")
            org_count = line[2]
            if org_count not in known_statistics:
                known_statistics[org_count] = 1
            else:
                known_statistics[org_count] += 1
        for i in known_statistics:
            known_names.append(i)
            known_counts.append(known_statistics[i])
        my_circle = plt.Circle((0, 0), 0.7, color='white')
        #plt.pie(known_counts, labels=known_names, wedgeprops={'linewidth': 7, 'edgecolor': 'white'}, colors=['red', 'blue', 'skyblue', 'green', 'black', 'gold', 'violet', 'teal', 'orange', 'olive'])
        plt.pie(known_counts, labels=known_names, wedgeprops={'linewidth': 7, 'edgecolor': 'white'}, colors=['red', 'blue', 'skyblue', 'green', 'black', 'gold', 'violet', 'teal', 'orange', 'olive'],normalize=True)
        plt.title("Shared Resistance Genes")
        p = plt.gcf()
        p.gca().add_artist(my_circle)
        plt.savefig(os.path.join(maindir, "res.jpg"))
        p.clf()
    if os.path.exists(bgc_path):
        for line in bgc_file:
            if line.startswith("#"):
                continue
            line = line.split("\t")
            org_count = len(line[1].split(","))
            if org_count not in bgc_statistics:
                bgc_statistics[org_count] = 1
            else:
                bgc_statistics[org_count] += 1
        for i in bgc_statistics:
            bgc_names.append(i)
            bgc_counts.append(bgc_statistics[i])
        my_circle = plt.Circle((0, 0), 0.7, color='white')
        #plt.pie(bgc_counts, labels=bgc_names, wedgeprops={'linewidth': 7, 'edgecolor': 'white'}, colors=['red', 'blue', 'skyblue', 'green', 'black', 'gold', 'violet', 'teal', 'orange', 'olive'])
        plt.pie(bgc_counts, labels=bgc_names, wedgeprops={'linewidth': 7, 'edgecolor': 'white'}, colors=['red', 'blue', 'skyblue', 'green', 'black', 'gold', 'violet', 'teal', 'orange', 'olive'],normalize=True)
        plt.title("Shared BGC Counts")
        p = plt.gcf()
        p.gca().add_artist(my_circle)
        plt.savefig(os.path.join(maindir, "bgc.jpg"))
        p.clf()
    if os.path.exists(core_path):
        for line in core_file:
            if line.startswith("#"):
                continue
            line = line.split("\t")
            org_count = len(line[8].split(","))
            if org_count not in core_statistics:
                core_statistics[org_count] = 1
            else:
                core_statistics[org_count] += 1
        for i in core_statistics:
            core_names.append(i)
            core_counts.append(core_statistics[i])
        my_circle = plt.Circle((0, 0), 0.7, color='white')
        plt.pie(core_counts, labels=core_names, wedgeprops={'linewidth': 7, 'edgecolor': 'white'}, colors=['red', 'blue', 'skyblue', 'green', 'black', 'gold', 'violet', 'teal', 'orange', 'olive'])
        plt.title("Shared Core Genes")
        p = plt.gcf()
        p.gca().add_artist(my_circle)
        plt.savefig(os.path.join(maindir, "core.jpg"))
        p.clf()
    if os.path.exists(dup_path):
        for line in dup_file:
            if line.startswith("#"):
                continue
            line = line.split("\t")
            org_count = line[2]
            if org_count not in dup_statistics:
                dup_statistics[org_count] = 1
            else:
                dup_statistics[org_count] += 1
        for i in dup_statistics:
            dup_names.append(i)
            dup_counts.append(dup_statistics[i])
        my_circle = plt.Circle((0, 0), 0.7, color='white')
        plt.pie(dup_counts, labels=dup_names, wedgeprops = { 'linewidth' : 7, 'edgecolor' : 'white' }, colors=['red', 'blue', 'skyblue', 'green', 'black', 'gold', 'violet', 'teal', 'orange', 'olive'])
        plt.title("Shared Duplication")
        p = plt.gcf()
        p.gca().add_artist(my_circle)
        plt.savefig(os.path.join(maindir,"duplication.jpg"))

def combine_bigscape_results(bigscape_result_dir, outdir, regions_to_filename_dict):
    # with open(os.path.join(bigscape_result_dir,"genbankdict.pickle"), "rb") as input_file:
    #     filename_dict = pickle.load(input_file)
    #Generate this from regions_to_filename_dict
    filename_dict = {os.path.splitext(os.path.basename(k))[0]:[v] for k,v in regions_to_filename_dict.items()}
    combined_log.info("FNdict: %s"%filename_dict)
    if os.path.exists(outdir) != True:
        os.mkdir(outdir)
    if os.path.exists(os.path.join(outdir, "tables")) != True:
        os.mkdir(os.path.join(outdir, "tables"))
    json_dict = {}
    json_result_file = os.path.join(os.path.join(outdir, "tables"), "combined_bgc_table.json")
    combined_table_file = open(os.path.join(os.path.join(outdir, "tables"),  "combined_bgc_table.tsv"), "w")
    network_files_folder = os.path.join(bigscape_result_dir, "network_files")
    if os.path.exists(network_files_folder) != True:
        combined_log.error("Bigscape result folder not found")
        return
    network_folders_list = os.listdir(network_files_folder)
    if len(network_folders_list) > 1:
        folder_name = ""
        combined_log.critical("Lethal warning: More than one network folder")
    else:
        folder_name = network_folders_list[0]
    bigscape_result_folder_path = os.path.join(network_files_folder, folder_name)
    network_annotations_file = open(os.path.join(bigscape_result_folder_path, "Network_Annotations_Full.tsv"), "r")
    line_num = 0
    network_dict = {}
    for line in network_annotations_file:
        line_num += 1
        if line_num == 1:
            continue
        else:
            split_line = line.split("\t")
            bgc = split_line[0]
            network_dict[bgc] = split_line[1:]
    mix_folder = os.path.join(bigscape_result_folder_path, "mix")
    for file in os.listdir(mix_folder):
        if file.startswith("mix_clustering"):
            cluster_filepath = os.path.join(mix_folder, file)
    cluster_file = open(cluster_filepath, "r")
    cluster_dict = {}
    line_num = 0
    for line in cluster_file:
        line_num += 1
        if line_num == 1:
            continue
        else:
            split_line = line.split("\t")
            bgc_name = split_line[0]
            fam_num = split_line[-1].strip("\n")
            if fam_num not in cluster_dict:
                cluster_dict[fam_num] = [bgc_name]
            else:
                cluster_dict[fam_num].append(bgc_name)
    identifier_line = "#Family\tRegion\tProduct\tClass\tDescription\tFilename\n"
    combined_table_file.write(identifier_line)
    for family in cluster_dict:
        line_to_write = ""
        line_to_write = line_to_write + family + "\t" + str(cluster_dict[family]) + "\t"
        family_products = []
        family_classes = []
        family_descs = []
        family_filenames = []
        for region in cluster_dict[family]:
            family_filenames.append(filename_dict[region][0])
        for region in cluster_dict[family]:
            family_products.append(network_dict[region][2])
            family_classes.append(network_dict[region][3])
            family_descs.append(network_dict[region][1])
        line_to_write = line_to_write + str(family_products) + "\t" + str(family_classes) + "\t" + str(family_descs) + "\t" + str(family_filenames) + "\n"
        combined_table_file.write(line_to_write)
    combined_table_file.close()
    with open(json_result_file, 'w') as outfile:
        json.dump(json_dict, outfile, indent = 2)


##for rows in bigscape analysis -to be continued-##
def parse_json(result_dict,outdir, regions_to_filename_dict):
    all_locus_to_region_dict = {}
    for res_dir in result_dict:
        with open(os.path.join(res_dir,"locus_to_region.pickle"), "rb") as input_file:
            locus_to_region_dict = pickle.load(input_file)
        for i in locus_to_region_dict:
            all_locus_to_region_dict[i] = locus_to_region_dict[i]
    #res_inv = {v: k for k, v in result_dict.iteritems()}
    res_inv = {v: k for k, v in result_dict.items()}
    filename_dict = {os.path.splitext(os.path.basename(k))[0]: [v] for k, v in regions_to_filename_dict.items()}
    #file_inv = {v[0]: k for k, v in filename_dict.iteritems()}
    file_inv = {v[0]: k for k, v in filename_dict.items()}
    # print all_locus_to_region_dict.keys()[0]
    # print all_locus_to_region_dict
    # print res_inv
    # print file_inv
    prox_res_coregene_to_hit_dict = {}
    combined_table_file = open(os.path.join(os.path.join(outdir, "tables"),  "combined_bgc_table.tsv"), "r")
    new_combined_table_file = open(os.path.join(os.path.join(outdir, "tables"),  "new_combined_bgc_table.tsv"), "w")
    reg_to_cluster = {}
    for line in combined_table_file:
        cores = []
        res = []
        reg_to_cluster = {}
        if line.startswith("#"):
            new_combined_table_file.write("#Family\tMembers\tRegion\tProduct\tClass\tDescription\tFilename\tShared_core_num\tShared_res_num \tCore_hits\tResistance_hits\tReg_Clust_num\n")
            continue
        line2 = line.split("\t")
        reg_list = ast.literal_eval(line2[1])
        file_list = ast.literal_eval(line2[-1])
        core_hit_list = []
        res_hit_list = []
        full_clust_list = []
        for i in range(0, len(reg_list)):
            core_2 = []
            res_2 = []
            region = reg_list[i]
            filename = file_list[i]
            base_name = region.split(".region")[0]
            clust =  int(region.split(".region")[1])
            specific_res_folder = res_inv[filename]
            specific_as_name = os.path.join(specific_res_folder, "antismash", base_name)
            if specific_as_name not in all_locus_to_region_dict:
                specific_as_name = specific_as_name[0:-2]
            reg_num = all_locus_to_region_dict[specific_as_name]
            full_clust_name = "cluster-"+ str(reg_num) + "_" + str(clust)
            full_clust_list.append(full_clust_name)
            core_json_path = os.path.join(specific_res_folder, "tables", "bgctable.json")
            if os.path.exists(core_json_path) == True:
                with open(core_json_path) as f:
                    bgc_json = json.load(f)
                if full_clust_name in bgc_json:
                    hits = bgc_json[full_clust_name]["hits"]
                    if len(hits) == 0:
                        core = "no_hit"
                        res = "no_hit"
                        print (core_json_path, full_clust_name, "NO HIT DETECTED")
                    else:
                        for j in hits:
                            if j[4] == "Core":
                                core_2.append(j[1])
                            if j[4] == "ResModel":
                                res_2.append(j[1])
            core_hit_list.append(core_2)
            res_hit_list.append(res_2)
        core_shared = set(core_hit_list[0])
        for k in core_hit_list[1:]:
            core_shared.intersection_update(k)
        res_shared = set(res_hit_list[0])
        for k in res_hit_list[1:]:
            res_shared.intersection_update(k)
        line_to_write = "\t".join([line2[0], str(len(reg_list)), line2[1], line2[2], line2[3], line2[4], line2[5].strip("\n") , str(len(core_shared)), str(len(res_shared)), str(core_hit_list), str(res_hit_list), str(full_clust_list)]) + "\n"
        new_combined_table_file.write(line_to_write)
    new_combined_table_file.close()









if __name__ == '__main__':
    parse_json({'/home/turgutmesut/test_36_24/test_36_24_0': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_1': '/home/turgutmesut/example_2.gbff'}
               ,'/home/turgutmesut/test_36_24',{'/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region011.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region016.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region006.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region005.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region008.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region012.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region009.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region003.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region018.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region004.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region014.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region017.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region001.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region007.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region010.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region015.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region002.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_0/antismash/CP000667.1.region013.gbk': '/home/turgutmesut/example.gbk', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region001.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region008.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region026.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region006.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region020.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region021.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region027.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region022.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region016.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region011.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region003.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region019.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region010.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region017.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region002.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region023.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region015.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region013.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region007.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region018.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region028.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region014.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region004.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region009.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region005.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region025.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region024.gbk': '/home/turgutmesut/example_2.gbff', '/home/turgutmesut/test_36_24/test_36_24_1/antismash/NC_009953.1.region012.gbk': '/home/turgutmesut/example_2.gbff'})
    print("arts output path and input dict")
    # generate_plots("/media/ecoli/1tb/bigscape-master-new/arts_integrate_test_2_/")
    #parse_json({'/media/ecoli/1tb/positive_examples/example/gcf_names/arts_results/arts_results_2': '/media/ecoli/1tb/positive_examples/example/gcf_names/GCF_000783435.2_ASM78343v2_genomic.gbff', '/media/ecoli/1tb/positive_examples/example/gcf_names/arts_results/arts_results_3': '/media/ecoli/1tb/positive_examples/example/gcf_names/GCF_000959485.1_ASM95948v1_genomic.gbff', '/media/ecoli/1tb/positive_examples/example/gcf_names/arts_results/arts_results_0': '/media/ecoli/1tb/positive_examples/example/gcf_names/GCF_000067205.1_ASM6720v1_genomic.gbff', '/media/ecoli/1tb/positive_examples/example/gcf_names/arts_results/arts_results_1': '/media/ecoli/1tb/positive_examples/example/gcf_names/GCF_000521505.1_ASM52150v1_genomic.gbff', '/media/ecoli/1tb/positive_examples/example/gcf_names/arts_results/arts_results_4': '/media/ecoli/1tb/positive_examples/example/gcf_names/GCF_001523725.1_ASM152372v1_genomic.gbff', '/media/ecoli/1tb/positive_examples/example/gcf_names/arts_results/arts_results_5': '/media/ecoli/1tb/positive_examples/example/gcf_names/GCF_900176645.1_LMG28154v1_genomic.gbff'})
