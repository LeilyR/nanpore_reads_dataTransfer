import rich_click as click
import os
import re
import sys
import yaml

def parsemodels(modelfile):
    dataStatus = False
    _lis = []
    with open(modelfile) as f:
        for line in f:
            if dataStatus:
                _split = re.sub(' +', ' ', line.strip()).split(' ')
                if len(_split) == 4:
                    _lis.append(
                        [
                            _split[0],
                            _split[1],
                            _split[2]
                        ]
                    )
                elif len(_split) == 5: # included barcodes.
                    _lis.append(
                        [
                            _split[0],
                            _split[1],
                            _split[3]
                        ]
                    )
            if line.strip().startswith('flowcell'):
                dataStatus = True
    return(_lis)

def modellist_to_dict(_models, modeldir):
    modeldic = {}
    for l in _models:
        flowcell, kit, model = l
        if flowcell not in modeldic:
           modeldic[flowcell] = {}
        # Assume every model has a high accuracy mode.
        if 'hac' in model:
            # for promethion runs, this would the priority list:
            # 1. model_sup_prom.cfg
            # 2. model_sup.cfg
            # 3. model_hac_prom.cfg
            if 'prom' in model:
                supmod_prom = os.path.join(
                    modeldir,
                    model.replace('hac', 'sup') + '.cfg'
                )
                supmod_reg = os.path.join(
                    modeldir,
                    model.replace('hac', 'sup').replace('_prom', '') + '.cfg'
                )
                if os.path.exists(supmod_prom):
                    supmod = supmod_prom
                elif os.path.exists(supmod_reg):
                    supmod = supmod_reg
            else:
                supmod = os.path.join(
                    modeldir,
                    model.replace('hac', 'sup') + '.cfg'
                )
            if os.path.exists(supmod):
                if kit in modeldic[flowcell]:
                    bps_present = int(modeldic[flowcell][kit].split('_')[3].replace('bps', ''))
                    bps_new = int(model.split('_')[3].replace('bps', ''))
                    if bps_new > bps_present:
                        print('replacing {} with {}'.format(modeldic[flowcell][kit], supmod))
                        modeldic[flowcell][kit] = supmod
                else:
                    modeldic[flowcell][kit] = supmod
            else:
                if kit in modeldic[flowcell]:
                    bps_present = int(modeldic[flowcell][kit].split('_')[3].replace('bps', ''))
                    bps_new = int(model.split('_')[3].replace('bps', ''))
                    if bps_new > bps_present:
                        print('replacing {} with {}'.format(modeldic[flowcell][kit], supmod))
                        modeldic[flowcell][kit] = os.path.join(modeldir, model +'.cfg')
                else:
                    modeldic[flowcell][kit] = os.path.join(modeldir, model +'.cfg')
        else:
            modeldic[flowcell][kit] = os.path.join(modeldir, model + '.cfg')
    return(modeldic)

@click.command(
    context_settings=dict(
        help_option_names=["-h", "--help"]
    )
)
@click.option(
   "-m",
   "--models",
   type=click.Path(exists=True),
   required=True,
   help='Specify a flowcell-kit-model file'
)
@click.option(
    '-d',
    '--modeldir',
    type=click.Path(exists=True),
    required=True,
    help='Specify the path to basecaller data directory (datadir, not the bindir).'
)
@click.option(
    '-o',
    '--outputdir',
    type=click.Path(exists=True),
    required=False,
    help='Specify the output directory to write ont_models.yaml into.',
    default=os.path.expanduser('~/configs/')
)
def main(models, modeldir,outputdir):
    if not os.path.exists(
        modeldir
    ):
        sys.exit('{} not found.'.format(modeldir))
    _models = parsemodels(models)
    modeldic = modellist_to_dict(_models, modeldir)
    ofile = os.path.join(outputdir, 'ont_models.yaml')
    with open(ofile, 'w') as f:
        yaml.dump(modeldic, f, default_flow_style=False)

