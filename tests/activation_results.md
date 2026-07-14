# Activation results (3 blind LLM routers, majority vote)
Generated 2026-07-13 from the build.py-generated trigger description (803 chars).
```
prompt                 exp  votes  fire? route  lang   component(maj)
------------------------------------------------------------------------------
markers_scrnaseq       True 3/3    ✓     ✓      ✓      ['cosg']
cosg_no_branding       True 3/3    ✓     ✓      ✓      ['cosg']
laris_named            True 3/3    ✓     ✓      ·      ['laris']
emergene_named         True 3/3    ✓     ✓      ·      ['emergene']
annotate_h5ad          True 3/3    ✓     ✓      ✓      ['piaso']
score_gene_set         True 3/3    ✓     ✓      ·      ['piaso']
lr_single_cell         True 3/3    ✓     ✓      ·      ['scalar']
lr_spatial_visium      True 3/3    ✓     ✓      ·      ['laris']
fast_rank_genes        True 3/3    ✓     ✓      ·      ['cosg']
cosg_seurat_R          True 3/3    ✓     ✓      ✓      ['cosgr']
cosg_h5ad_py           True 3/3    ✓     ✓      ✓      ['cosg']
neg_nfcore             False 0/3    ✓     ·      ·      [None]
neg_qc                 False 0/3    ✓     ·      ·      [None]
neg_scvi               False 0/3    ✓     ·      ·      [None]
------------------------------------------------------------------------------
Activation (fire/no-fire): 14/14
Component routing (when firing): 11/11
Language routing (probes): 5/5
RESULT: ALL PASS
```

Notes: one router individually routed 'score a gene set' to Emergene rather than PIASO
(both expose a gene-set scoring function); majority was PIASO. All 3 routers unanimously
declined the 3 negative controls (nf-core, plain QC, scVI batch correction).
