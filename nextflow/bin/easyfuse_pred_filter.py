#!/usr/bin/env python

import pandas as pd
import sys

easyfuse_pred = sys.argv[1]
easyfuse_pred_filter = sys.argv[2]

easyfuse_pred_pd = pd.read_csv(easyfuse_pred, sep=';', low_memory=False)
easyfuse_pred_pd['ft_junc_span_best'] = easyfuse_pred_pd['ft_junc_best'] + easyfuse_pred_pd['ft_span_best']

easyfuse_pred_pd_groupby = easyfuse_pred_pd[['Fusion_Gene','ft_junc_span_best']].groupby(by = 'Fusion_Gene', as_index = False).max()
easyfuse_pred_pd_read = pd.merge(easyfuse_pred_pd,easyfuse_pred_pd_groupby,on=['Fusion_Gene','ft_junc_span_best'],how='inner').drop_duplicates()
easyfuse_pred_pd_read_groupby = easyfuse_pred_pd_read[['Fusion_Gene','prediction_prob']].groupby(by = 'Fusion_Gene', as_index = False).max()
easyfuse_pred_pd_read_prob = pd.merge(easyfuse_pred_pd_read,easyfuse_pred_pd_read_groupby,on=['Fusion_Gene','prediction_prob'],how='inner').drop_duplicates()

easyfuse_pred_pd_read_prob_frame = easyfuse_pred_pd_read_prob[easyfuse_pred_pd_read_prob['frame']!='no_frame']

easyfuse_pred_pd_read_prob_frame.to_csv(easyfuse_pred_filter, header=True, index=False, sep='\t')
