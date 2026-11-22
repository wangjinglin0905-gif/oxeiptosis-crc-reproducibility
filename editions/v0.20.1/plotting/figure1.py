"""Build the AJCR composite figure from the Valdeolivas section ZIP and OX4 v0.20.1.

Usage: python figure1.py SECTION_ZIP UNPACKED_OX4_0_20_1 OUTPUT_DIR [ARIAL_FONT_DIR]
Requires Python, NumPy, Pillow and Arial regular/bold. Does not require network access.
"""
import csv
import gzip
import hashlib
import io
import json
import math
import sys
import zipfile
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

SAMPLE = 'SN048_A121573_Rep1'
GENES = ('KEAP1', 'PGAM5', 'AIFM1', 'OTUD1')

def digest(file):
    return hashlib.sha256(file.read_bytes()).hexdigest()

def read_points(archive, release):
    frozen = release/'historical_spatial/outputs/Zenodo7760264_gate_v01_2026-11-17'
    with (frozen/'barcode_label_coordinate_join.csv').open(encoding='utf-8-sig') as handle:
        rows = {row['barcode']:row for row in csv.DictReader(handle) if row['sample']==SAMPLE}
    with (release/'tables/S33_section_results.csv').open(encoding='utf-8-sig') as f:
        sections = {r['sample']:r for r in csv.DictReader(f) if r['construct']=='OX4'}
    assert SAMPLE==min(sections)
    with zipfile.ZipFile(archive) as z:
        def member(rel): return z.read(SAMPLE+'/'+rel)
        def matrix_text(file):
            return gzip.decompress(member('filtered_feature_bc_matrix/'+file+'.gz')).decode().splitlines()
        png=member('spatial/tissue_hires_image.png')
        raw=Image.open(io.BytesIO(png)).convert('RGB')
        scale=json.loads(member('spatial/scalefactors_json.json'))
        pos={r[0]:r for r in csv.reader(member('spatial/tissue_positions_list.csv').decode().splitlines())}
        barcodes=matrix_text('barcodes.tsv')
        features=[l.split('\t') for l in matrix_text('features.tsv')]
        assert {r[2] for r in features}=={'Gene Expression'}
        target={}
        for i,g in enumerate(GENES):
            ids=[j+1 for j,r in enumerate(features) if r[1]==g]
            assert len(ids)==1
            target[ids[0]]=i
        counts=np.zeros((len(barcodes),4),dtype=np.int64)
        totals=np.zeros(len(barcodes),dtype=np.int64)
        with z.open(SAMPLE+'/filtered_feature_bc_matrix/matrix.mtx.gz') as file, gzip.open(file,'rb') as stream:
            header=next(line for line in stream if not line.startswith(b'%'))
            ng,nb,nnz=map(int,header.split())
            assert (ng,nb)==(len(features),len(barcodes))
            seen=0
            for line in stream:
                gene,spot,value=map(int,line.split())
                totals[spot-1]+=value
                if gene in target:counts[spot-1,target[gene]]+=value
                seen+=1
            assert seen==nnz
    keep=[];points=[]
    for j,b in enumerate(barcodes):
        r=rows[b]; orig=pos[b]
        assert int(r['total_gene_expression_counts'])==totals[j] and int(r['in_tissue'])==int(orig[1])
        if orig[1]!='1' or totals[j]<500 or r['coarse_label'] in ('author_exclude','unlabelled'):continue
        x=float(orig[5])*scale['tissue_hires_scalef']; y=float(orig[4])*scale['tissue_hires_scalef']
        assert abs(float(r['x_hires'])-x)<.01 and abs(float(r['y_hires'])-y)<.01
        assert 0<=x<raw.width and 0<=y<raw.height
        keep.append(j)
        points.append([b,r['coarse_label'],x,y])
    v=np.log1p(counts[keep,:]*10000/totals[keep,None])
    scores=((v-v.mean(0))/v.std(0,ddof=1)).mean(1)
    scores/=scores.std(ddof=1)
    for row,score in zip(points,scores):row.append(float(score))
    t=np.array([p[-1] for p in points if p[1]=='tumor'])
    s=np.array([p[-1] for p in points if p[1]=='stroma'])
    delta=float(np.sign(t[:,None]-s[None,:]).mean())
    assert (len(points),len(t),len(s))==(2049,337,249)
    assert abs(delta-float(sections[SAMPLE]['cliff_delta']))<1e-12
    info=dict(sample=SAMPLE,image_sha256=hashlib.sha256(png).hexdigest(),eligible_spots=len(points),
              tumour=len(t),stroma=len(s),cliff_delta=delta,score_above_2=int(sum(scores>2)),
              score_below_minus2=int(sum(scores< -2)),hires_scalefactors=scale)
    return raw,points,info

def colour(v):
    stops=[(43,40,112),(52,109,151),(41,151,132),(142,193,85),(250,222,56)]
    f=min(4,max(0,v+2));j=min(3,int(f))
    return tuple(round(a+(b-a)*(f-j)) for a,b in zip(stops[j],stops[j+1]))

def render(archive,release,out,fontdir):
    out.mkdir(parents=True,exist_ok=True)
    im,pts,info=read_points(archive,release)
    with (out/'Figure_1_spots.csv').open('w',encoding='utf8',newline='') as f:
        w=csv.writer(f);w.writerow(['barcode','source_compartment_group','x_hires','y_hires','OX4']);w.writerows(pts)
    baseplot=release/'figure_sources/01_quantitative_panels_v0.20.0.tif'
    assert digest(baseplot)=='5ba7a988b74ddf3a75397bb4c87ee1389904ed45399db4a3e39357e935d137ff'
    quant=Image.open(baseplot).convert('RGB')
    assert quant.size==(5860,1899)
    font=lambda px,bold=False:ImageFont.truetype(str(fontdir/('arialbd.ttf' if bold else 'arial.ttf')),px)
    lower=ImageDraw.Draw(quant)
    for label,x in zip('def',[63,1969,3933]):
        lower.rectangle([x-4,12,x+70,180],fill='white')
        lower.text((x-4,49),label,font=font(101,True),fill='black',anchor='lt')
    width=3900
    quant=quant.resize((width,1264),Image.Resampling.LANCZOS)
    canvas=Image.new('RGB',(width,2865),'white')
    canvas.paste(quant,(0,1599))
    p=ImageDraw.Draw(canvas)
    scale=info['hires_scalefactors']
    diam=scale['spot_diameter_fullres']*scale['tissue_hires_scalef']
    x,y=np.array([pt[2:4] for pt in pts]).T
    side=math.ceil(max(np.ptp(x),np.ptp(y))+3*diam)
    x0=math.floor((min(x)+max(x)-side)/2);y0=math.floor((min(y)+max(y)-side)/2)
    base=im.crop((x0,y0,x0+side,y0+side))
    size=1190
    tones={'tumor':(0,101,151),'stroma':(210,94,15)}
    starts=[57,1358,2659]
    for k,(left,title) in enumerate(zip(starts,['H&E','Histological compartments','OX4 transcript score'])):
        p.text((left,28),chr(97+k),font=font(66,True),fill='black')
        p.text((left+88,34),title,font=font(58),fill='black')
        art=base.copy() if k==0 else Image.blend(base,Image.new('RGB',base.size,'white'),.57)
        draw=ImageDraw.Draw(art)
        if k:
            for b,label,cx,cy,v in pts:
                if k==1 and label not in tones:continue
                cx-=x0;cy-=y0
                draw.ellipse((cx-diam/2,cy-diam/2,cx+diam/2,cy+diam/2),fill=tones[label] if k==1 else colour(v))
        canvas.paste(art.resize((size,size),Image.Resampling.LANCZOS),(left,139))
    # 55 micrometre Visium capture spots; section-specific Space Ranger calibration.
    bar=round(500/55*diam*size/side)
    bar_rects=[]
    for left in starts:
        xbar=left+13
        rect=(xbar,1350,xbar+bar-1,1362)
        p.rectangle(rect,fill='black')
        p.text((xbar+bar+20,1325),'500 μm',font=font(48),fill='black')
        bar_rects.append(rect)
    for j,(lab,c) in enumerate(tones.items()):
        lx=starts[1]+j*485
        p.ellipse([lx,1452,lx+36,1488],fill=c)
        p.text((lx+51,1432),('Tumour','Stroma')[j],font=font(50),fill='black')
    barx=starts[2]+28
    for dx in range(670):
        p.line((barx+dx,1450,barx+dx,1483),fill=colour(-2+4*dx/669))
    for tick,xp in zip(['≤−2','0','≥2'],[barx,barx+335,barx+670]):
        p.text((xp,1491),tick,font=font(47),anchor='mt',fill='black')
    p.text((barx+750,1445),'Section SD',font=font(48),fill='black')
    canvas.save(out/'Figure_1.tif',compression='tiff_lzw',dpi=(600,600))
    info.update(archive_sha256=digest(archive),source_doi='10.5281/zenodo.7760264',source_license='CC BY 4.0',
                crop_hires_pixels=[x0,y0,side,side],scale_bar_um=500,spot_diameter_um=55,
                scale_bar_width_export_pixels=bar,scale_bar_rects_export_pixels=bar_rects,
                lower_panels_source_sha256=digest(baseplot),output_tiff_sha256=digest(out/'Figure_1.tif'),
                source_selection='First alphanumeric eligible section ID; no OX4 values used for selection.')
    (out/'Figure_1_provenance.json').write_text(json.dumps(info,indent=2),encoding='utf8')
    print('Figure 1',len(pts),info['cliff_delta'],flush=True)

if __name__=='__main__':
    render(Path(sys.argv[1]),Path(sys.argv[2]),Path(sys.argv[3]),Path(sys.argv[4]) if len(sys.argv)>4 else Path('C:/Windows/Fonts'))
