#include "graph.h"
#include "utility/graph_operations.h"

Graph::Graph(){
   	
    vertices_count_ = 0;
    edges_count_ = 0;
    labels_count_ = 0;
    max_degree_ = 0;
    max_label_frequency_ = 0;

    offsets_ = NULL;
    neighbors_ = NULL;
    labels_ = NULL;

    reverse_index_offsets_ = NULL;
    reverse_index_ = NULL;
    labels_frequency_.clear();
#if OPTIMIZED_LABELED_GRAPH == 1
    nlf_ = NULL;
#endif
}

Graph::~Graph() {
    delete[] offsets_;
    delete[] neighbors_;
    delete[] labels_;
    delete[] reverse_index_offsets_;
    delete[] reverse_index_;
#if OPTIMIZED_LABELED_GRAPH == 1
    delete[] nlf_;
#endif
}

void Graph::loadGraphFromFile(const std::string &file_path) {
    std::ifstream infile(file_path);

    if (!infile.is_open()) {
        std::cout << "Can not open the graph file " << file_path << " ." << std::endl;
        exit(-1);
    }

    char type;
    infile >> type >> vertices_count_ >> edges_count_;
    offsets_ = new ui[vertices_count_ + 1];
    offsets_[0] = 0;

    neighbors_ = new VertexID[edges_count_ * 2];
    labels_ = new LabelID[vertices_count_];

    LabelID max_label_id = 0;
    std::vector<ui> neighbors_offset(vertices_count_, 0);

    while (infile >> type) {
        if (type == 'v') { // Read vertex.
            VertexID id;
            LabelID  label;
            ui degree;
            infile >> id >> label >> degree;

            labels_[id] = label;
            offsets_[id + 1] = offsets_[id] + degree;

            if (degree > max_degree_) {
                max_degree_ = degree;
            }

            if (labels_frequency_.find(label) == labels_frequency_.end()) {
                labels_frequency_[label] = 0;
                if (label > max_label_id)
                    max_label_id = label;
            }

            labels_frequency_[label] += 1;
        }
        else if (type == 'e') { // Read edge.
            VertexID begin;
            VertexID end;
            infile >> begin >> end;

            ui offset = offsets_[begin] + neighbors_offset[begin];
            neighbors_[offset] = end;

            offset = offsets_[end] + neighbors_offset[end];
            neighbors_[offset] = begin;

            neighbors_offset[begin] += 1;
            neighbors_offset[end] += 1;
        }
    }

    infile.close();
    labels_count_ = (ui)labels_frequency_.size() > (max_label_id + 1) ? (ui)labels_frequency_.size() : max_label_id + 1;

    for (auto element : labels_frequency_) {
        if (element.second > max_label_frequency_) {
            max_label_frequency_ = element.second;
        }
    }

    for (ui i = 0; i < vertices_count_; ++i) {
        std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1]);
    }

    BuildReverseIndex();

#if OPTIMIZED_LABELED_GRAPH == 1
    BuildNLF();
#endif

}

void Graph::loadGraphFromFileCompressed(const std::string &degree_path, const std::string &edge_path,
                                        const std::string &label_path) {
    
}

void Graph::printGraphMetaData() {
    std::cout << "|V|: " << vertices_count_ << ", |E|: " << edges_count_ << ", |\u03A3|: " << labels_count_ << std::endl;
    std::cout << "Max Degree: " << max_degree_ << ", Max Label Frequency: " << max_label_frequency_ << std::endl;
}

void Graph::buildCoreTable() {
    core_table_ = new int[vertices_count_];
    GraphOperations::getKCore(this, core_table_);

    for (ui i = 0; i < vertices_count_; ++i) {
    	// 2-core
        if (core_table_[i] > 1) {
            core_length_ += 1;
        }
    }
}

void Graph::BuildReverseIndex() {
    reverse_index_ = new ui[vertices_count_];
    reverse_index_offsets_= new ui[labels_count_ + 1];
    reverse_index_offsets_[0] = 0;

    ui total = 0;
    // reverse_index_offsets_[0] = 0;
    // reverse_index_offsets_[1] = 0;
    // reverse_index_offsets_[i + 2] = reverse_index_offsets_[i + 1] + labels_frequency_[i];
    // reverse_index_offsets_[i] stores the start position of vertex with label i - 1;
    for (ui i = 0; i < labels_count_; ++i) {
        reverse_index_offsets_[i + 1] = total;
        total += labels_frequency_[i];
    }

    // write vertices and update reverse_index_offsets_ array so that reverse_index_offsets_[i] stores the start position of vertex with label i;
    for (ui i = 0; i < vertices_count_; ++i) {
        LabelID label = labels_[i];
        reverse_index_[reverse_index_offsets_[label + 1]++] = i;
    }
}

#if OPTIMIZED_LABELED_GRAPH == 1
void Graph::BuildNLF() {
    nlf_ = new std::unordered_map<LabelID, ui>[vertices_count_];
    for (ui i = 0; i < vertices_count_; ++i) {
        ui count;
        const VertexID * neighbors = getVertexNeighbors(i, count);

        for (ui j = 0; j < count; ++j) {
            VertexID u = neighbors[j];
            LabelID label = getVertexLabel(u);
            if (nlf_[i].find(label) == nlf_[i].end()) {
                nlf_[i][label] = 0;
            }

            nlf_[i][label] += 1;
        }
    }
}
#endif