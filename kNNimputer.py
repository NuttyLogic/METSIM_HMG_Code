#! /usr/env python3

import numpy as np
import pandas as pd
import random


class kNNimputer:

    def __init__(self, dataframe=None, distance='euclidean'):
        self.df = dataframe
        self.distance = distance
        self.row_labels = self.df.index
        self.samples = list(self.df)
        self.values = np.transpose(self.df.values)
        self.distance_matrix = None
        self.chromosome_pairwise_dict = None
        self.neighbor_dicts = None

    def run(self, imputation_dist=6000000, boundary_dist=2000000, k=5,
            missing_value_tolerance=0.0, variance_tolerance=None):
        print('Getting Global Neighbors')
        global_neighbors = kNNimputer.euclidean(self.values)
        print('Processing Windows')
        self.chromosome_pairwise_windows(imputation_dist=imputation_dist,
                                         boundary_dist=boundary_dist,
                                         global_neighbors=global_neighbors)
        print('Getting ' + str(k) + ' Nearest Neighbors')
        self.get_nearest_neighbors(k=k)
        print('Imputing Missing Values')
        self.nearest_neighbor_imputation(missing_value_tolerance=missing_value_tolerance,
                                         variance_tolerance=variance_tolerance)

    @staticmethod
    def euclidean(array=None, global_neighbors=None):
        """ Pairwise euclidean distance
        ----------------------------------------
        input: numpy array with samples as rows
        output: pairwise euclidean distance array"""
        # split numpy array by columns
        methylation_samples = np.split(array, 1)
        # output to ensure correct matrix orientation
        # initialize empty list to hold pairwise distance values
        pairwise_distance = []
        # iterate through sample combinations
        for count, methylation_vector in enumerate(methylation_samples[0]):
            # add list for all sample pairwise comparison
            pairwise_distance.append([])
            for comparison_count, comparison_vector in enumerate(methylation_samples[0]):
                # only perform calculation for upper half of distance matrix to save time
                if comparison_count >= count:
                    # stack sample values
                    comparison_array = np.stack((methylation_vector, comparison_vector), axis=0)
                    # drop any comparison value with NA for any sample, note this is only done pairwise between
                    # samples to get best distance comparison estimate
                    comparison_array = np.ma.compress_cols(np.ma.masked_invalid(comparison_array))
                    try:
                        pairwise_distance[count].append(np.linalg.norm(comparison_array[0] - comparison_array[1]))
                    # if all values are missing the missing neighbor is close,
                    except IndexError:
                        if global_neighbors:
                            pairwise_distance[count].append(global_neighbors[count][comparison_count])
                        else:
                            pairwise_distance[count].append(0)
                else:
                    pairwise_distance[count].append(0)
        # add upper and lower matrices
        pairwise_upper = np.asarray(pairwise_distance)
        pairwise_lower = np.transpose(pairwise_upper)
        # return list to ensure compatibility with downstream tools
        pairwise_distance = (pairwise_upper + pairwise_lower).tolist()
        return pairwise_distance

    @staticmethod
    def chromosome_imputation_chunks(site_labels, imputation_distance=3000000, boundary=1000000):
        """Segments chromosome inputs into discreet chunks for imputation. Three overlapping windows traverse the
        genome. Sites are assigned to imputation window if they reside in the middle section of the relevant imputation
        window. Sites at beginning or end of a chromosome are assigned to the first or last window accordingly even
        though the site is outside the window inner boundary.
        --------------------------
        input; site_labels, list of chromosome labels, chr1 or 1, with : separating genome coordinates, ie 1:100000
        input; imputation_distance=3000000, size of sliding window
        input; boundary=1000000, size of inner boundaries in sliding window, ie a 1,000,000 boundary with a 3,000,000bp
            sliding window will result in a window segmented into three 1,000,000bp chunks
        returns; window_list; a list of sites to use for pairwise distance calculation
        returns; window_dict; dict for every site listing the correct window to use for imputation"""
        # initialize placeholder values
        placeholder_list = [None, None, None, None, None]
        windows_list = []
        windows = [[], [], []]
        window_count = 0
        site_window_dict = {}
        # iterate through site labels
        for site in site_labels:
            # split site by position and chromosome
            site_split = site.split(':')
            site_chromosome = site_split[0]
            site_position = int(site_split[1])
            # initialize variable for first site
            if not placeholder_list[0]:
                placeholder_list = [site_chromosome, site_position, site_position + boundary,
                                    site_position + 2 * boundary, site_position + imputation_distance]
                windows[0].append(site)
                # assign site to imputation window
                site_window_dict[site] = window_count
            # set rules for chromosome change
            elif placeholder_list[0] != site_chromosome:
                placeholder_list = [site_chromosome, site_position, site_position + boundary,
                                    site_position + 2 * boundary, site_position + imputation_distance]
                windows_list.append(windows[0])
                windows_list.append(windows[0])
                windows_list.append(windows[0])
                windows = [[site], [], []]
                # advance window count
                window_count += 1
                # set dictionary value for new window
                site_window_dict[site] = window_count
            elif placeholder_list[0] == site_chromosome:
                # rules to add sites to current window
                # if site location not past second internal boundary add to current imputation window
                if site_position <= placeholder_list[2]:
                    windows[0].append(site)
                    # rule to capture early sites
                    site_window_dict[site] = window_count
                elif placeholder_list[2] < site_position <= placeholder_list[3]:
                    # middle sites are in original window and next window
                    windows[0].append(site)
                    windows[1].append(site)
                    site_window_dict[site] = window_count
                elif placeholder_list[3] < site_position <= placeholder_list[4]:
                    # if site is past the second internal boundary site imputed at next window
                    windows[0].append(site)
                    windows[1].append(site)
                    windows[2].append(site)
                    site_window_dict[site] = window_count + 1
                elif site_position > placeholder_list[4]:
                    # initialize next window
                    windows_list.append(windows[0])
                    windows = [list(windows[1]), list(windows[2]), []]
                    window_count += 1
                    placeholder_list = [site_chromosome, placeholder_list[2], placeholder_list[3],
                                        placeholder_list[3] + boundary,
                                        placeholder_list[1] + boundary + imputation_distance]
                    # place out of bounds site in current windows
                    if site_position <= placeholder_list[2]:
                        windows[0].append(site)
                        site_window_dict[site] = window_count
                    elif placeholder_list[2] < site_position <= placeholder_list[3]:
                        windows[0].append(site)
                        windows[1].append(site)
                        site_window_dict[site] = window_count
                    elif placeholder_list[3] < site_position <= placeholder_list[4]:
                        windows[0].append(site)
                        windows[1].append(site)
                        windows[2].append(site)
                        site_window_dict[site] = window_count + 1
                    elif site_position > placeholder_list[4]:
                        placeholder_list = [site_chromosome, site_position, site_position + boundary,
                                            site_position + 2 * boundary, site_position + imputation_distance]
                        windows_list.append(windows[0])
                        windows_list.append(windows[0])
                        windows_list.append(windows[0])
                        windows = [[site], [], []]
                        # advance window count
                        window_count += 1
                        # set dictionary value for new window
                        site_window_dict[site] = window_count
        # after loop add all windows to list for imputation
        windows_list.append(windows[0])
        windows_list.append(windows[0])
        windows_list.append(windows[0])
        return windows_list, site_window_dict

    def chromosome_pairwise_windows(self, imputation_dist=6000000, boundary_dist=2000000, global_neighbors=None):
        """Calculate sliding imputation windows. Each site is assigned to the imputation window where the site lies on
         in the middle section
         -----------------------------
         input; self.row_labels, list of genome coordinate seperated by :, ie. chr1:11111111
         returns; self.distance_matrix, list of sites to include in pairwise distance matrix
         returns; self.chromosome_pairwise_dict, hash table linking coordinates to imputation matrix"""
        windows_list, site_window_dict = kNNimputer.chromosome_imputation_chunks(self.row_labels,
                                                                                 imputation_distance=imputation_dist,
                                                                                 boundary=boundary_dist)
        pairwise_arrays = []
        for window in windows_list:
            window_df = self.df.loc[window]
            window_values = np.transpose(window_df.values)
            pairwise_array = kNNimputer.euclidean(window_values, global_neighbors)
            pairwise_arrays.append(pairwise_array)
        self.distance_matrix = pairwise_arrays
        self.chromosome_pairwise_dict = site_window_dict

    def get_nearest_neighbors(self, k=5):
        neighbor_dicts = []
        for distance_matrix in self.distance_matrix:
            local_neighbor_dict = {}
            for sample in zip(self.samples, distance_matrix):
                sample_key = list(self.samples)
                sorted_sample_key = [a for b, a in sorted(zip(sample[1], sample_key))]
                neighbors = [sample_key.index(x) for x in sorted_sample_key[1:(k+1)]]
                local_neighbor_dict[sample_key.index(sample[0])] = neighbors
            neighbor_dicts.append(local_neighbor_dict)
        self.neighbor_dicts = neighbor_dicts

    def pairwise_matrix(self):
        self.distance_matrix = [kNNimputer.euclidean(self.values)]

    def nearest_neighbor_imputation(self, missing_value_tolerance=0.9, variance_tolerance=None):
        # drop low information rows
        numeric_count = self.df.count(axis=1) / len(list(self.df))
        # Pandas series object, drop categories below missing_value_tolerance
        below_threshold_index = numeric_count.drop(numeric_count[numeric_count >= missing_value_tolerance].index)
        # Designate columns above variance_drop_proportion threshold
        # of series labels
        self.df = self.df.drop(below_threshold_index.index)
        # if the columns are labeled by a number pandas will pull by the label instead of the index which can lead to
        # unpredictable results, ensure all labels are strings to get around this
        self.df.columns = [str(x) for x in self.samples]
        for index, row in self.df.iterrows():
            null_methylation_values = np.argwhere(np.isnan(row.values))
            if len(null_methylation_values) > 0:
                if self.chromosome_pairwise_dict:
                    null_site_neighbor_dict = self.neighbor_dicts[self.chromosome_pairwise_dict[index]]
                else:
                    null_site_neighbor_dict = self.neighbor_dicts[0]
                for null_site in null_methylation_values:
                    neighbors = null_site_neighbor_dict[null_site[0]]
                    if variance_tolerance:
                        if (row[neighbors]).var(skipna=True) > variance_tolerance:
                            print('what')
                    neighbor_average = round((row[neighbors]).mean(skipna=True), 3)
                    row[null_site[0]] = neighbor_average


def accuracy_assessment(test_df=None, proportion_of_random_sites=0.01, k=5,
                        imputation_dist=6000000, boundary_dist=2000000):
    # drop all rows with missing information
    test_df = test_df[test_df.notnull().all(1)]
    validation_df = pd.DataFrame.copy(test_df)
    ix = [(row, col) for row in range(test_df.shape[0]) for col in range(test_df.shape[1])]
    for row, col in random.sample(ix, int(round(proportion_of_random_sites * len(ix)))):
        test_df.iat[row, col] = np.nan
    nan_index = np.where(np.asanyarray(np.isnan(test_df)))
    test_imputer = kNNimputer(dataframe=test_df)
    test_imputer.run(imputation_dist=imputation_dist, boundary_dist=boundary_dist, k=k, missing_value_tolerance=0.0)
    # absolute_error, validation_value, test_value, cpg_site_index
    validation_list = [[], [], [], []]
    for index in zip(nan_index[0], nan_index[1]):
        validation_list[3].append(index[0])
        test_value = test_imputer.df.iat[index[0], index[1]]
        validation_value = validation_df.iat[index[0], index[1]]
        absolute_error = np.absolute(test_value - validation_value)
        validation_list[2].append(test_value)
        validation_list[0].append(absolute_error)
        validation_list[1].append(validation_value)
    return validation_list
