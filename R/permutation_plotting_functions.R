get_topology = function(graph) {
    adjacency_matrix = convert_to_adjacency_matrix(graph)
    topological_table = matrix(rep(0, 7), nrow = 1, ncol = 7)
    colnames(topological_table) = c("Number of vertices", "Number of edges", "Density", "Centralization", "Heterogeneity", "Mean clustering coefficient",
        "Mean scaled connectivity")

    # Estimation of the network topological parameters
    statistics = fundamentalNetworkConcepts(adjacency_matrix, GS = NULL)
    topological_table[, 1] = vcount(graph)
    topological_table[, 2] = ecount(graph)
    topological_table[, 3] = statistics$Density
    topological_table[, 4] = statistics$Centralization
    topological_table[, 5] = statistics$Heterogeneity
    topological_table[, 6] = mean(statistics$ClusterCoef)
    topological_table[, 7] = mean(statistics$ScaledConnectivity)
    return(topological_table)
}

seed_and_extend = function(algorithm, seeds, disease_name, no_of_permutations = 1, hub_names, network, type, number_of_hubs,
    network_vertices, network_weights, output_directory) {
    results_from_subnetworks = list(hub_frequency_table = matrix(0, ncol = length(extract_percent), nrow = number_of_hubs, dimnames = list(hub_names,
        paste(extract_percent, "%", sep = ""))), performance_with_seeds = NULL, performance_without_seeds = NULL)
    observ_pred_df = set_observed_pred_table(network_vertices, seeds)

    for (percent in extract_percent) {
        for (i in 1:no_of_permutations) {
            # Sublist extraction
            percent_sublist = unique(sample(seeds, round(length(seeds) * percent/100)))

            # Check if the substructed list has minimal size (>1)
            if (length(percent_sublist) > 1) {
                constructed_network = construct_network(algorithm, network, percent_sublist, network_weights)
                results_from_subnetworks = evaluate_networks(constructed_network, hub_names, seeds, network_vertices, percent,
                  observ_pred_df, percent_sublist, i, results_from_subnetworks)
            } else {
                break
            }
        }
        results_from_subnetworks[[1]][, paste(percent, "%", sep = "")] = results_from_subnetworks[[1]][, paste(percent, "%",
            sep = "")]/no_of_permutations * 100
    }
    write.table(results_from_subnetworks[[1]], file = file.path(output_directory, paste0(disease_name, "_", type, "_", no_of_permutations,
        "_permut_hub_frequency_", algorithm, ".txt")), sep = "\t", quote = FALSE)
    write.table(results_from_subnetworks[[2]], file = file.path(output_directory, paste0(disease_name, "_", type, "_", no_of_permutations,
        "_permut_performance_with_seeds_", algorithm, ".txt")), sep = "\t", quote = FALSE)
    write.table(results_from_subnetworks[[3]], file = file.path(output_directory, paste0(disease_name, "_", type, "_", no_of_permutations,
        "_permut_performance_without_seeds_", algorithm, ".txt")), sep = "\t", quote = FALSE)
}

evaluate_networks = function(constructed_network, hub_names, seeds, network_vertices, percent, observ_pred_df, percent_sublist,
    i, results_from_subnetworks) {
    hub_frequency_table = hub_frequency(hub_names, percent, constructed_network, results_from_subnetworks[[1]])
    subnetwork_vertices = vertex.attributes(constructed_network)[[1]]

    # Create permutation dataframe with true and predicted values
    list_observ_pred = true_vs_predicted_df(network_vertices, seeds, subnetwork_vertices, percent_sublist, observ_pred_df)

    # Create four-table with respective calculations
    current_permut = paste(percent, "%", " (permut ", i, ")", sep = "")
    performance_with_seeds = calculate_performance(list_observ_pred[[1]], results_from_subnetworks[[2]])
    row.names(performance_with_seeds)[nrow(performance_with_seeds)] = current_permut

    performance_without_seeds = calculate_performance(list_observ_pred[[2]], results_from_subnetworks[[3]])
    row.names(performance_without_seeds)[nrow(performance_without_seeds)] = current_permut
    return(list(hub_frequency_table, performance_with_seeds, performance_without_seeds))
}

true_vs_predicted_df = function(network_vertices, seeds, subnetwork_vertices, percent_sublist, observ_pred_df) {
    observ_pred_df[which(rownames(observ_pred_df) %in% subnetwork_vertices), 2] = 1

    observ_pred_df_removed = observ_pred_df
    for (gene in percent_sublist) {
        observ_pred_df_removed = observ_pred_df_removed[-which(rownames(observ_pred_df_removed) == gene), ]
    }
    return(list(observ_pred_df, observ_pred_df_removed))
}

construct_network = function(algorithm, graph, percent_sublist, network_weights) {
  constructed_network = NULL

  if (length(grep("Markov_walk", algorithm)) == 1) {
    steps = strsplit(algorithm, "_")[[1]][1]
    constructed_network = Markov_walk(graph, percent_sublist, steps)
  } else {
    constructed_network = switch(algorithm,
                                 first_neighborhood = first_neighborhood(graph, percent_sublist),
                                 shortest_path = shortest_path(graph, percent_sublist, network_weights))
  }
  return(constructed_network)
}

calculate_performance = function(observ_pred_df, performance) {
    ## predAct is two col dataframe of pred,act
    trues = observ_pred_df[, 1]
    preds = observ_pred_df[, 2]
    xTab = table(preds, trues)
    if (nrow(xTab) == 1) {
        xTab = rbind(xTab, c(0, 0))
    }
    class = as.character(sort(unique(preds)))
    add = matrix(NA, ncol = 5, nrow = 1, dimnames = list(c(), c("FP rate", "TP rate(Recall)", "Precision", "Accuracy", "F-score")))

    add[1, 1] = xTab[2, 1]/sum(xTab[, 1])  #FP rate
    add[1, 2] = xTab[2, 2]/sum(xTab[, 2])  #TP rate
    add[1, 3] = xTab[2, 2]/sum(xTab[2, ])  #Precision
    add[1, 4] = sum(xTab[2, 2], xTab[1, 1])/sum(xTab)  #Accuracy
    add[1, 5] = add[1, 3] * add[1, 2]  #F-score

    performance = rbind(performance, add)
    return(performance)
}

hub_frequency = function(hub_names, percent, constructed_network, hub_frequency_table) {
    for (hub in hub_names) {
        if (hub %in% V(constructed_network)$name) {
            hub_frequency_table[hub, paste(percent, "%", sep = "")] = hub_frequency_table[hub, paste(percent, "%", sep = "")] +
                1
        }
    }
    return(hub_frequency_table)
}

process_resulting_files = function(file_name, extract_percent, no_of_permutations = 1) {
    performance = read.table(file_name, sep = "\t", row.names = 1)
    summary_table = NULL
    median_performance = NULL
    n = 1
    for (i in 1:length(extract_percent)) {
        sds = apply(performance[n:(n + no_of_permutations - 1), ], 2, sd)
        med = colMedians(performance[n:(n + no_of_permutations - 1), ])
        if (!(is.na(med[1]))) {
            summary_table = rbind(summary_table, paste(round(med, 2), round(sds, 2), sep = " +/- "))
            row.names(summary_table)[nrow(summary_table)] = strsplit(rownames(performance[n, ]), " ")[[1]][1]

            median_performance = rbind(median_performance, med)
            row.names(median_performance)[nrow(median_performance)] = paste(strsplit(rownames(performance[n, ]), " ")[[1]][1],
                " median", sep = "")
        }
        n = n + no_of_permutations
    }
    AUC = get_AUC(median_performance)
    return(list(median_performance, AUC, summary_table))
}

get_AUC = function(median_performance) {
    FP = c(0, median_performance[, 1], 1)
    TP = c(0, median_performance[, 2], 1)

    # AUC = sum((FP[2:length(TP)]-FP[1:length(TP)-1])*TP[2:length(TP)])
    AUC = round(as.double((FP[2:length(FP)] - FP[2:length(FP) - 1]) %*% (TP[2:length(FP)] + TP[2:length(FP) - 1]))/2, digits = 2)
    return(AUC)
}

plot_degree_frequency = function(extract_percent, hub_frequency_tables, output_directory, disease_name, algorithms) {
    number_of_hubs = 5
    colors = matrix(brewer.pal(length(algorithms), "Set1"), ncol = 1)
    rownames(colors) = algorithms

    if (length(grep("degree_discounted", names(hub_frequency_tables))) >= 1) {
        png(file = file.path(output_directory, paste0(disease_name, "_hubs_discounted.png")), height = 5, width = 8.5, units = "in",
            res = 300)
        par(mar = c(5.1, 4.1, 4.1, 17), xpd = TRUE)
    } else {
        png(file = file.path(output_directory, paste0(disease_name, "_hubs.png")), height = 5, width = 7, units = "in", res = 300)
        par(mar = c(5.1, 4.1, 4.1, 13), xpd = TRUE)
    }
    for (i in 1:length(hub_frequency_tables)) {
        algorithm = unlist(strsplit(names(hub_frequency_tables)[i], " "))[2]
        line_type = 1
        color = colors[algorithm, 1]
        if (length(grep("degree_discounted", names(hub_frequency_tables)[i])) >= 1) {
            line_type = 3
        }
        if (!(is.null(hub_frequency_tables[[i]]))) {
            mean_frequency = colMeans(hub_frequency_tables[[i]][1:number_of_hubs, ])
            if (i == 1) {
                plot(extract_percent, mean_frequency, type = "l", lty = line_type, col = color, ylim = c(0, 100), ylab = "Mean frequency of top 5 hubs (%)",
                  xlab = "Percent of input seeds")
            } else {
                lines(extract_percent, mean_frequency, lty = line_type, col = color)
            }
        }
    }
    legend("right", inset = c(-0.7, 0), algorithms, col = colors[, 1], bty = "n", pch = 19, title = "Algorithms", y.intersp = 2)
    dev.off()
}

plot_curves = function(median_performance, extract_percent, disease_name, algorithms, output_directory, name) {
    AUC_legend = NULL
    colors = matrix(brewer.pal(length(algorithms), "Set1"), ncol = 1)
    rownames(colors) = algorithms

    # Plotting into separate png files
    if (length(grep("degree_discounted", names(median_performance))) > 0) {
        file_name = file.path(output_directory, paste0(disease_name, "_ROC_", name, "_discounted.png"))
    } else {
        file_name = file.path(output_directory, paste0(disease_name, "_ROC_", name, ".png"))
    }
    png(file = file_name, height = 7, width = 8.5, units = "in", res = 300)
    par(mar = c(5.1, 4.1, 4.1, 13), xpd = TRUE)

    for (i in 1:length(median_performance)) {
        # TP rate and FP rate curve
        FP = median_performance[[i]][[1]][, 1]
        TP = median_performance[[i]][[1]][, 2]

        algorithm = unlist(strsplit(names(median_performance)[i], " "))[2]
        AUC_legend = c(AUC_legend, paste(algorithm, "=", median_performance[[i]][[2]]))

        line_type = 1
        shape = 20
        color = colors[algorithm, 1]
        if (length(grep("degree_discounted", names(median_performance)[i])) >= 1) {
            line_type = 3
            shape = 4
        }
        if (i == 1) {
            plot(c(0, FP, 1), c(0, TP, 1), col = color, type = "l", lty = line_type, xlab = "FP rate", ylab = "TP rate", main = paste(disease_name,
                name))
        } else {
            lines(c(0, FP, 1), c(0, TP, 1), lty = line_type, col = color)
        }
        points(FP, TP, col = color, pch = shape)
    }
    legend("right", inset = c(-0.5, 0), AUC_legend, col = colors, bty = "n", pch = 19, title = "AUCs of algorithms", y.intersp = 2)
    dev.off()

    legend = NULL
    if (length(grep("degree_discounted", names(median_performance))) >= 1) {
        file_name = file.path(output_directory, paste0(disease_name, "_prec-recall_", name, "_discounted.png"))
    } else {
        file_name = file.path(output_directory, paste0(disease_name, "_prec-recall_", name, ".png"))
    }
    png(file = file_name, height = 7, width = 8.5, units = "in", res = 300)
    par(mar = c(5.1, 4.1, 4.1, 13), xpd = TRUE)

    for (i in 1:length(median_performance)) {
        # Precision and recall curve
        recall = median_performance[[i]][[1]][, 2]
        precision = median_performance[[i]][[1]][, 3]
        algorithm = unlist(strsplit(names(median_performance)[i], " "))[2]
        legend = c(legend, algorithm)

        line_type = 1
        shape = 20
        color = colors[algorithm, 1]
        if (length(grep("degree_discounted", names(median_performance)[i])) >= 1) {
            line_type = 3
            shape = 4
        }
        if (i == 1) {
            plot(recall, precision, col = color, type = "l", lty = line_type, xlab = "recall", ylab = "precision", xlim = c(0,
                1), ylim = c(0, 1), main = paste(disease_name, name))
        } else {
            lines(recall, precision, col = color, lty = line_type)
        }
        points(recall, precision, col = color, pch = shape)
    }
    legend("right", inset = c(-0.4, 0), legend, col = colors, bty = "n", pch = 19, title = "Algorithms", y.intersp = 2)
    dev.off()
}

plot_AUCs = function(AUCs_with_seeds, AUCs_without_seeds, output_directory, algorithms, type) {
    colors = brewer.pal(ncol(AUCs_with_seeds), "Set1")
    shapes = c(15, 16, 17, 7, 8, 4, 18)

    png(file = file.path(output_directory, paste0("AUCs_", type, ".png")), height = 5, width = 6.5, units = "in", res = 300)
    par(mar = c(5.1, 4.1, 4.1, 13), xpd = TRUE)
    for (i in 1:nrow(AUCs_with_seeds)) {
        if (i == 1) {
            plot(AUCs_with_seeds[i, ], AUCs_without_seeds[i, ], xlim = c(0.4, 1), ylim = c(0.4, 1), pch = shapes[i], xlab = "AUC seeds included",
                ylab = "AUC seeds excluded", col = colors)
        } else {
            points(AUCs_with_seeds[i, ], AUCs_without_seeds[i, ], pch = shapes[i], col = colors)
        }
    }
    legend("bottomright", inset = c(-0.69, -0.2), colnames(AUCs_with_seeds), col = "black", bty = "n", pch = shapes, title = "Algorithm shapes",
        y.intersp = 1.2)
    legend("topright", inset = c(-0.8, -0.2), rownames(AUCs_with_seeds), col = colors, bty = "n", pch = 19, title = "Disease colors",
        y.intersp = 1.2)
    dev.off()
}
