#! /usr/env python3

import numpy as np
import pandas as pd
import scipy.stats as stats
from sklearn import preprocessing


class DataCleanup:
    """Note! This class is designed for use with phenotype data, to normalize information for feature selection. Whole,
    data table will be loaded into memory so this is not suitable for extremely large data sets"""

    def __init__(self, file, separator, header=True, named_rows=True):
        """Simple dataframe import for normalization and data cleaning
        ---------------------------------------------
        file = file to import, csv or tsv
        separator = file feature used to designate columns
        header = True: if true header will be imported as the first row
        named_rows = True: first column represents row names
        missing_value_representation = how missing values are represented in the data table, default is blank"""
        if header and named_rows:
            imported_data_frame = pd.DataFrame.from_csv(file, sep=separator, header=0, index_col=0)
        elif header and not named_rows:
            imported_data_frame = pd.DataFrame.from_csv(file, sep=separator, header=0)
        elif not header and named_rows:
            imported_data_frame = pd.DataFrame.from_csv(file, sep=separator, index_col=0)
        else:
            imported_data_frame = pd.DataFrame.from_csv(file, sep=separator)
        self.df = imported_data_frame
        self.dropped_row_index = None
        self.dropped_categories = None

    @staticmethod
    def rank_based_inv_norm(x, c=0):
        """
        Perform the rank-based inverse normal transformation on data x.

        Reference: Beasley TM, Erickson S, Allison DB. Rank-based inverse normal
        transformations are increasingly used, but are they merited? Behav Genet.
        2009; 39(5):580-95.

        :param x: input array-like
        :param c: constant parameter
        :return: y
        """
        return stats.norm.ppf((stats.rankdata(x) + c) / (x.size - 2 * c + 1))

    def list_column_removal(self, removal_list=None):
        """Removes columns supplied in list
        ---------------------------------------------
        removal_list=None: list of columns to remove
        """
        self.df = self.df.drop(removal_list, axis=1)

    def list_row_removal(self, removal_list=None):
        """Removes rows supplied in list
        ---------------------------------------------
        removal_list=None: list of rows to remove from the dataframe
        """
        self.df = self.df.drop(removal_list, axis=0)

    def numeric_data_transformation(self, errors='coerce'):
        """Transform variables to a numeric representation, any non-numeric values will be represented as
        np.nan (not a number).
        ---------------------------------------------
        errors='coerce'; how to handle missing or non-numeric data in pandas; options
            ‘raise’, then invalid parsing will raise an exception
            ‘coerce’, then invalid parsing will be set as NaN
            ‘ignore’, then invalid parsing will return the input"""
        # transform value to numeric representation, if not a numeric value will be replaced with np.nan
        self.df = self.df.apply(pd.to_numeric, errors=errors)

    def remove_low_information_rows(self, missing_value_tolerance=0.5):
        """Remove rows with a lot of missing information, important for large data sets that require imputation
        ---------------------------------------------
        missing_value_tolerance=0.5, 0-1 value: rows where the proportion of missing values is below
            this value will be removed from the dataframe"""
        # Proportion of column with values missing
        numeric_count = self.df.count(axis=1) / len(list(self.df))
        # Pandas series object, drop categories below missing_value_tolerance
        below_threshold_index = numeric_count.drop(numeric_count[numeric_count >= missing_value_tolerance].index)
        # Designate columns above variance_drop_proportion threshold, filtered_categories.axes[0].tolist() returns list
        # of series labels
        self.df = self.df.drop(list(below_threshold_index.index))
        self.dropped_row_index = list(below_threshold_index.index)

    def remove_low_information_columns(self, missing_value_tolerance=0.8):
        """Transform variables to a numeric representation, any non-numeric values will be represented as
        np.nan (not a number).
        ---------------------------------------------
        missing_value_tolerance=0.8, 0-1 value: categories where the proportion of numeric values is below
        this value will be removed from the dataframe, this includes all non-numeric columns"""
        # Proportion of column with values missing
        numeric_count = self.df.count(axis=0) / len(self.df.index)
        # Pandas series object, drop categories below missing_value_tolerance
        above_threshold_categories = numeric_count.drop(numeric_count[numeric_count <= missing_value_tolerance].index)
        # Designate columns above variance_drop_proportion threshold, filtered_categories.axes[0].tolist() returns list
        # of series labels
        self.df = self.df[above_threshold_categories.axes[0].tolist()]
        self.dropped_categories = above_threshold_categories

    def standard_imputation(self, imputation_strategy='most_frequent'):
        """Imputation using mean, meadian, or most-frequent values. Standard sklearn implementation wraped around a
        pandas data frame
        ---------------------------------------------
        imputation_strategy='most_frequent', sklearn allows for 'mean', 'median', or the 'most_frequent' value to
            be used for imputation, choose the necessary imputation technique"""
        # set sklearn imputation object, and set axis=0 to work with column formatted data
        data_frame_imputer = preprocessing.Imputer(missing_values=np.nan, strategy=imputation_strategy, axis=0)
        # Fit the imputation model and transform data
        imputed_data_frame = data_frame_imputer.fit_transform(self.df)
        # Re-initialize dataframe with imputed data
        self.df = pd.DataFrame(imputed_data_frame, index=list(self.df.index), columns=list(self.df))

    def normal_transformation(self, normal_test=1000, normal_percentile=0.01):
        """Perform a test of normality, if the test is below threshold perform rank based inverse transformation. If
        above threshold, standardize data with a mean of 0 and stdev of 1.
        ---------------------------------------------
        normal_tests=1000, sets the probability cutoff normal distributions by generating random samples of the length
            of the input data and scoring them for normality, default is 1000 random samples
        normal_percentile=0.1, percentile in the distribution to used to set your normal distribution threshold before
            transformation"""
        # generate randomly sampled data from normal distribution and score for normality
        # storage list
        test_scores = []
        # generate random distribution the same length as the input data
        test_length = len(list(self.df.index))
        for _ in range(normal_test):
            # generate random distributions with a mean of 0 and stdev of 1
            random_normal = np.random.normal(loc=0, scale=1, size=test_length)
            chi_squared_stat, two_sided_chi_p_value = stats.normaltest(random_normal)
            test_scores.append(two_sided_chi_p_value)
        test_scores.sort()
        # set normal threshold
        percentile_cutoff_index = int(round(normal_percentile * len(test_scores)))
        normal_threshold = test_scores[percentile_cutoff_index]
        for column in list(self.df):
            # preform scipy.stats.normaltest
            chi_squared_stat, two_sided_chi_p_value = stats.normaltest(self.df[column])
            # if distribution score below threshold rank based inverse normal transformation else scale (mean=0, std=1)
            if two_sided_chi_p_value <= normal_threshold:
                self.df[column] = DataCleanup.rank_based_inv_norm(self.df[column])
            else:
                self.df[column] = preprocessing.scale(self.df[column])

    def remove_low_variance_features(self, variance_threshold=0.0):
        """This function should only be called on normalized and imputed dataframe, drop columns with
        low variance that are likely non-informative.
        ---------------------------------------------
        variance_threshold=0.0; threshold to drop columns
        """
        variance = self.df.var(axis=0)
        above_threshold_categories = variance.drop(variance[variance <= variance_threshold].index)
        self.df = self.df[above_threshold_categories.index.tolist()]

    def export_normalized_data(self, output_path=None, separator=None):
        """Export dataframe
        ---------------------------------------------
        output_path=None: path to export file
        separator=None: whitespace indicator"""
        self.df.to_csv(path_or_buf=output_path, sep=separator, header=True, index=True)
