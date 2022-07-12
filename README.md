# visual\_fq

A fastq file visualization tool 🦀📈🚀

## Dependencies

### Ubuntu Linux:

```sudo apt install pkg-config libfreetype6-dev libfontconfig1-dev```

### CentOS Linux:

```sudo yum install pkg-config freetype-devel fontconfig-devel```

## install

install rust first

```bash
git clone https://github.com/sharkLoc/visual_fq.git
cd visual_fq
cargo b --release
```

## usage

```bash
./target/release/visual_fq

[Usage information]
uncompressed fastq:
        ./target/release/visual_fq your_file.fq > out.matrix 2> summary.log

compressed fastq:
         gzip -dc your_file.fq.gz | ./target/release/visual_fq /dev/stdin > out.matrix 2> summary.log

result files:
         out.matrix  summary.log  Base_plot.png
```

## output
output plot: `Base_plot.png`

![Base_plot](https://user-images.githubusercontent.com/50580507/176699721-8045127f-6dab-4caf-98a7-fb62f758fd3d.png)

output plot: `Qual_plot.png`

![Qual_plot](https://user-images.githubusercontent.com/50580507/178418104-c11a2a7c-59cc-47eb-9a85-fa085db8a1ec.png)


output file: `out.matrix` 
```
Iterm	Total	A	T	G	C	N	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41
pos:1	100000	30.70	25.84	21.82	21.63	0.01	0	0	13	0	0	0	0	0	0	0	0	0	1430	0	0	0	0	0	0	0	0	0	0	0	0	0	0	2624	0	0	0	0	95933	0	0	0	0	0	0	0	0	0
pos:2	100000	34.55	31.93	18.97	14.54	0.00	0	0	0	0	0	0	0	0	0	0	0	0	385	0	0	0	0	0	0	0	0	0	0	0	0	0	0	911	0	0	0	0	98703	0	0	0	0	1	0	0	0	0
pos:3	100000	32.02	33.43	17.12	17.44	0.00	0	0	0	0	0	0	0	0	0	0	0	0	877	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1946	0	0	0	0	15710	0	0	0	0	81467	0	0	0	0
pos:4	100000	31.65	31.75	17.97	18.63	0.00	0	0	0	0	0	0	0	0	0	0	0	0	780	0	0	0	0	0	0	0	0	0	129	0	0	0	0	1596	0	0	0	0	6439	0	0	0	0	91056	0	0	0	0
pos:5	100000	33.47	30.47	18.11	17.96	0.00	0	0	0	0	0	0	0	0	0	0	0	0	17577	0	0	0	0	0	0	0	0	0	5466	0	0	0	0	6954	0	0	0	0	17572	0	0	0	0	52431	0	0	0	0
pos:6	100000	32.12	31.54	17.47	18.87	0.00	0	0	0	0	0	0	0	0	0	0	0	0	679	0	0	0	0	0	0	0	0	0	146	0	0	0	0	1277	0	0	0	0	6271	0	0	0	0	26316	0	0	0	65311
pos:7	100000	32.92	31.16	18.04	17.87	0.00	0	0	0	0	0	0	0	0	0	0	0	0	9002	0	0	0	0	0	0	0	0	0	2271	0	0	0	0	6266	0	0	0	0	12529	0	0	0	0	22133	0	0	0	47799
pos:8	100000	32.66	31.96	17.56	17.82	0.00	0	0	0	0	0	0	0	0	0	0	0	0	9214	0	0	0	0	0	0	0	0	0	1448	0	0	0	0	6513	0	0	0	0	13309	0	0	0	0	21591	0	0	0	47925
pos:9	100000	30.80	33.43	17.45	18.32	0.00	0	0	0	0	0	0	0	0	0	0	0	0	453	0	0	0	0	0	0	0	0	0	621	0	0	0	0	517	0	0	0	0	4449	0	0	0	0	17064	0	0	0	76896
pos:10	100000	31.46	32.12	18.02	18.40	0.00	0	0	0	0	0	0	0	0	0	0	0	0	491	0	0	0	0	0	0	0	0	0	575	0	0	0	0	787	0	0	0	0	2828	0	0	0	0	10635	0	0	0	84684
pos:11	100000	31.94	31.52	18.13	18.41	0.00	0	0	0	0	0	0	0	0	0	0	0	0	927	0	0	0	0	0	0	0	0	0	428	0	0	0	0	1398	0	0	0	0	3826	0	0	0	0	11294	0	0	0	82127
pos:12	100000	31.97	31.99	18.10	17.94	0.00	0	0	0	0	0	0	0	0	0	0	0	0	345	0	0	0	0	0	0	0	0	0	352	0	0	0	0	557	0	0	0	0	1985	0	0	0	0	8106	0	0	0	88655
pos:13	100000	32.51	31.52	17.98	17.99	0.00	0	0	0	0	0	0	0	0	0	0	0	0	695	0	0	0	0	0	0	0	0	0	759	0	0	0	0	1165	0	0	0	0	3008	0	0	0	0	10280	0	0	0	84093
pos:14	100000	32.29	31.67	18.02	18.02	0.00	0	0	0	0	0	0	0	0	0	0	0	0	607	0	0	0	0	0	0	0	0	0	656	0	0	0	0	1027	0	0	0	0	2771	0	0	0	0	9584	0	0	0	85355
pos:15	100000	32.52	31.59	18.25	17.64	0.00	0	0	0	0	0	0	0	0	0	0	0	0	5121	0	0	0	0	0	0	0	0	0	3160	0	0	0	0	4679	0	0	0	0	7625	0	0	0	0	17783	0	0	0	61632
pos:16	100000	32.29	31.70	18.00	18.01	0.00	0	0	0	0	0	0	0	0	0	0	0	0	432	0	0	0	0	0	0	0	0	0	496	0	0	0	0	575	0	0	0	0	3163	0	0	0	0	11843	0	0	0	83491
pos:17	100000	32.08	31.87	17.86	18.19	0.00	0	0	0	0	0	0	0	0	0	0	0	0	1083	0	0	0	0	0	0	0	0	0	1066	0	0	0	0	1407	0	0	0	0	4047	0	0	0	0	11861	0	0	0	80536
pos:18	100000	31.94	31.49	18.22	18.35	0.00	0	0	0	0	0	0	0	0	0	0	0	0	565	0	0	0	0	0	0	0	0	0	598	0	0	0	0	926	0	0	0	0	3010	0	0	0	0	10255	0	0	0	84646
pos:19	100000	32.82	30.86	19.40	16.92	0.00	0	0	0	0	0	0	0	0	0	0	0	0	12832	0	0	0	0	0	0	0	0	0	5492	0	0	0	0	7460	0	0	0	0	10008	0	0	0	0	19335	0	0	0	44873
pos:20	100000	32.66	31.23	18.47	17.64	0.00	0	0	0	0	0	0	0	0	0	0	0	0	9836	0	0	0	0	0	0	0	0	0	5133	0	0	0	0	4498	0	0	0	0	12902	0	0	0	0	21850	0	0	0	45781
pos:21	100000	32.09	31.68	18.13	18.10	0.00	0	0	0	0	0	0	0	0	0	0	0	0	894	0	0	0	0	0	0	0	0	0	1260	0	0	0	0	3062	0	0	0	0	3379	0	0	0	0	19436	0	0	0	71969
pos:22	100000	32.45	30.92	18.38	18.25	0.00	0	0	0	0	0	0	0	0	0	0	0	0	5324	0	0	0	0	0	0	0	0	0	3456	0	0	0	0	5914	0	0	0	0	6264	0	0	0	0	18111	0	0	0	60931
pos:23	100000	32.01	31.90	18.12	17.96	0.00	0	0	0	0	0	0	0	0	0	0	0	0	958	0	0	0	0	0	0	0	0	0	1233	0	0	0	0	2806	0	0	0	0	3264	0	0	0	0	15490	0	0	0	76249
pos:24	100000	31.86	31.93	18.00	18.22	0.00	0	0	0	0	0	0	0	0	0	0	0	0	3334	0	0	0	0	0	0	0	0	0	2580	0	0	0	0	4910	0	0	0	0	5495	0	0	0	0	17609	0	0	0	66072
pos:25	100000	32.16	31.51	18.16	18.17	0.00	0	0	0	0	0	0	0	0	0	0	0	0	3386	0	0	0	0	0	0	0	0	0	2702	0	0	0	0	4698	0	0	0	0	5190	0	0	0	0	17410	0	0	0	66614
pos:26	100000	31.75	31.51	18.68	18.05	0.00	0	0	0	0	0	0	0	0	0	0	0	0	10560	0	0	0	0	0	0	0	0	0	3764	0	0	0	0	6248	0	0	0	0	6106	0	0	0	0	17460	0	0	0	55862
pos:27	100000	31.84	31.70	18.14	18.32	0.00	0	0	0	0	0	0	0	0	0	0	0	0	1605	0	0	0	0	0	0	0	0	0	1599	0	0	0	0	3798	0	0	0	0	3168	0	0	0	0	16208	0	0	0	73622
pos:28	100000	31.90	31.67	18.25	18.18	0.00	0	0	0	0	0	0	0	0	0	0	0	0	2800	0	0	0	0	0	0	0	0	0	1657	0	0	0	0	3165	0	0	0	0	3459	0	0	0	0	14050	0	0	0	74869
pos:29	100000	32.14	31.45	19.21	17.20	0.00	0	0	0	0	0	0	0	0	0	0	0	0	20265	0	0	0	0	0	0	0	0	0	5506	0	0	0	0	8284	0	0	0	0	7202	0	0	0	0	17728	0	0	0	41015
pos:30	100000	31.68	31.89	18.26	18.16	0.00	0	0	0	0	0	0	0	0	0	0	0	0	5889	0	0	0	0	0	0	0	0	0	3956	0	0	0	0	7350	0	0	0	0	5600	0	0	0	0	21968	0	0	0	55237
pos:31	100000	31.89	31.18	18.40	18.53	0.00	0	0	0	0	0	0	0	0	0	0	0	0	2366	0	0	0	0	0	0	0	0	0	2110	0	0	0	0	4397	0	0	0	0	3951	0	0	0	0	17784	0	0	0	69392
pos:32	100000	31.91	31.16	18.56	18.37	0.00	0	0	0	0	0	0	0	0	0	0	0	0	1168	0	0	0	0	0	0	0	0	0	1037	0	0	0	0	2259	0	0	0	0	2713	0	0	0	0	13156	0	0	0	79667
pos:33	100000	31.66	31.72	18.74	17.88	0.00	0	0	0	0	0	0	0	0	0	0	0	0	11484	0	0	0	0	0	0	0	0	0	3726	0	0	0	0	6220	0	0	0	0	6213	0	0	0	0	17312	0	0	0	55045
pos:34	100000	31.56	31.63	18.59	18.22	0.00	0	0	0	0	0	0	0	0	0	0	0	0	2616	0	0	0	0	0	0	0	0	0	1987	0	0	0	0	4279	0	0	0	0	3711	0	0	0	0	16955	0	0	0	70452
pos:35	100000	31.97	31.56	18.03	18.43	0.00	0	0	0	0	0	0	0	0	0	0	0	0	673	0	0	0	0	0	0	0	0	0	710	0	0	0	0	1756	0	0	0	0	2038	0	0	0	0	11353	0	0	0	83470
pos:36	100000	31.91	31.47	18.24	18.38	0.00	0	0	0	0	0	0	0	0	0	0	0	0	3716	0	0	0	0	0	0	0	0	0	1808	0	0	0	0	3278	0	0	0	0	3826	0	0	0	0	13683	0	0	0	73689
pos:37	100000	31.88	30.89	18.61	18.62	0.00	0	0	0	0	0	0	0	0	0	0	0	0	8004	0	0	0	0	0	0	0	0	0	3065	0	0	0	0	5066	0	0	0	0	5129	0	0	0	0	15565	0	0	0	63171
pos:38	100000	31.91	31.64	19.00	17.44	0.00	0	0	0	0	0	0	0	0	0	0	0	0	21872	0	0	0	0	0	0	0	0	0	6000	0	0	0	0	8879	0	0	0	0	7359	0	0	0	0	18147	0	0	0	37743
pos:39	100000	31.72	31.66	18.15	18.46	0.00	0	0	0	0	0	0	0	0	0	0	0	0	1248	0	0	0	0	0	0	0	0	0	1801	0	0	0	0	4775	0	0	0	0	3383	0	0	0	0	22136	0	0	0	66657
pos:40	100000	31.70	31.79	18.30	18.22	0.00	0	0	0	0	0	0	0	0	0	0	0	0	1419	0	0	0	0	0	0	0	0	0	1094	0	0	0	0	2404	0	0	0	0	3151	0	0	0	0	14291	0	0	0	77641
pos:41	100000	31.49	31.76	18.32	18.44	0.00	0	0	0	0	0	0	0	0	0	0	0	0	2349	0	0	0	0	0	0	0	0	0	1383	0	0	0	0	2569	0	0	0	0	3247	0	0	0	0	12768	0	0	0	77684
pos:42	100000	31.35	31.67	18.50	18.48	0.00	0	0	0	0	0	0	0	0	0	0	0	0	1344	0	0	0	0	0	0	0	0	0	1028	0	0	0	0	2088	0	0	0	0	2496	0	0	0	0	10924	0	0	0	82120
pos:43	100000	31.50	31.67	18.60	18.22	0.00	0	0	0	0	0	0	0	0	0	0	0	0	1074	0	0	0	0	0	0	0	0	0	809	0	0	0	0	1655	0	0	0	0	2121	0	0	0	0	9821	0	0	0	84520
pos:44	100000	31.85	31.29	18.61	18.24	0.00	0	0	0	0	0	0	0	0	0	0	0	0	6144	0	0	0	0	0	0	0	0	0	2465	0	0	0	0	4173	0	0	0	0	4760	0	0	0	0	14513	0	0	0	67945
pos:45	100000	31.64	31.68	18.34	18.34	0.00	0	0	0	0	0	0	0	0	0	0	0	0	1287	0	0	0	0	0	0	0	0	0	1172	0	0	0	0	2378	0	0	0	0	2584	0	0	0	0	12642	0	0	0	79937
pos:46	100000	31.92	31.39	18.27	18.41	0.00	0	0	0	0	0	0	0	0	0	0	0	0	3820	0	0	0	0	0	0	0	0	0	2068	0	0	0	0	3682	0	0	0	0	4213	0	0	0	0	14494	0	0	0	71723
pos:47	100000	31.66	31.72	18.40	18.21	0.00	0	0	0	0	0	0	0	0	0	0	0	0	2190	0	0	0	0	0	0	0	0	0	1458	0	0	0	0	2864	0	0	0	0	3168	0	0	0	0	13478	0	0	0	76842
pos:48	100000	31.60	31.52	18.37	18.51	0.00	0	0	0	0	0	0	0	0	0	0	0	0	2498	0	0	0	0	0	0	0	0	0	1572	0	0	0	0	3014	0	0	0	0	3467	0	0	0	0	13210	0	0	0	76239
pos:49	100000	31.57	31.66	18.25	18.53	0.00	0	0	0	0	0	0	0	0	0	0	0	0	1915	0	0	0	0	0	0	0	0	0	1267	0	0	0	0	2447	0	0	0	0	2757	0	0	0	0	11881	0	0	0	79733
pos:50	100000	33.39	31.71	18.68	16.21	0.00	0	0	0	0	0	0	0	0	0	0	0	0	27146	0	0	0	0	0	0	0	0	0	6870	0	0	0	0	9544	0	0	0	0	7727	0	0	0	0	17074	0	0	0	31639
pos:51	100000	31.69	30.80	19.16	18.35	0.00	0	0	0	0	0	0	0	0	0	0	0	0	13613	0	0	0	0	0	0	0	0	0	5901	0	0	0	0	9789	0	0	0	0	7159	0	0	0	0	22990	0	0	0	40548
pos:52	100000	31.81	31.05	18.64	18.50	0.00	0	0	0	0	0	0	0	0	0	0	0	0	10511	0	0	0	0	0	0	0	0	0	5285	0	0	0	0	9518	0	0	0	0	6896	0	0	0	0	22259	0	0	0	45531
pos:53	100000	31.98	31.48	18.20	18.35	0.00	0	0	0	0	0	0	0	0	0	0	0	0	1714	0	0	0	0	0	0	0	0	0	1768	0	0	0	0	4414	0	0	0	0	4070	0	0	0	0	19177	0	0	0	68857
pos:54	100000	31.79	31.77	18.26	18.17	0.00	0	0	0	0	0	0	0	0	0	0	0	0	2742	0	0	0	0	0	0	0	0	0	1691	0	0	0	0	3308	0	0	0	0	4530	0	0	0	0	14971	0	0	0	72758
pos:55	100000	31.75	31.57	18.73	17.94	0.00	0	0	0	0	0	0	0	0	0	0	0	0	16997	0	0	0	0	0	0	0	0	0	4936	0	0	0	0	7804	0	0	0	0	7903	0	0	0	0	18106	0	0	0	44254
pos:56	100000	31.68	31.42	18.41	18.49	0.00	0	0	0	0	0	0	0	0	0	0	0	0	1555	0	0	0	0	0	0	0	0	0	1660	0	0	0	0	3859	0	0	0	0	4099	0	0	0	0	18966	0	0	0	69861
pos:57	100000	31.95	31.43	18.31	18.31	0.00	0	0	0	0	0	0	0	0	0	0	0	0	3980	0	0	0	0	0	0	0	0	0	2100	0	0	0	0	3953	0	0	0	0	5209	0	0	0	0	15842	0	0	0	68916
pos:58	100000	31.79	31.20	18.49	18.52	0.00	0	0	0	0	0	0	0	0	0	0	0	0	6545	0	0	0	0	0	0	0	0	0	2843	0	0	0	0	4981	0	0	0	0	5734	0	0	0	0	15826	0	0	0	64071
pos:59	100000	31.83	31.34	18.54	18.29	0.00	0	0	0	0	0	0	0	0	0	0	0	0	3089	0	0	0	0	0	0	0	0	0	2112	0	0	0	0	3925	0	0	0	0	4301	0	0	0	0	14809	0	0	0	71764
pos:60	100000	31.88	31.33	18.47	18.32	0.00	0	0	0	0	0	0	0	0	0	0	0	0	11941	0	0	0	0	0	0	0	0	0	4075	0	0	0	0	6582	0	0	0	0	6769	0	0	0	0	17330	0	0	0	53303
pos:61	100000	31.79	31.28	18.60	18.32	0.00	0	0	0	0	0	0	0	0	0	0	0	0	15541	0	0	0	0	0	0	0	0	0	5576	0	0	0	0	8508	0	0	0	0	8045	0	0	0	0	18901	0	0	0	43429
pos:62	100000	31.89	31.40	18.78	17.93	0.00	0	0	0	0	0	0	0	0	0	0	0	0	21292	0	0	0	0	0	0	0	0	0	6550	0	0	0	0	9736	0	0	0	0	8172	0	0	0	0	18198	0	0	0	36052
pos:63	100000	31.51	31.37	19.21	17.91	0.00	0	0	0	0	0	0	0	0	0	0	0	0	21733	0	0	0	0	0	0	0	0	0	7080	0	0	0	0	10933	0	0	0	0	8047	0	0	0	0	18576	0	0	0	33631
pos:64	100000	31.81	31.52	18.49	18.18	0.00	0	0	0	0	0	0	0	0	0	0	0	0	12778	0	0	0	0	0	0	0	0	0	6600	0	0	0	0	11152	0	0	0	0	7948	0	0	0	0	20927	0	0	0	40595
pos:65	100000	31.93	31.23	18.41	18.43	0.00	0	0	0	0	0	0	0	0	0	0	0	0	7266	0	0	0	0	0	0	0	0	0	4596	0	0	0	0	8283	0	0	0	0	10733	0	0	0	0	17673	0	0	0	51449
pos:66	100000	31.57	31.28	18.55	18.60	0.00	0	0	0	0	0	0	0	0	0	0	0	0	6755	0	0	0	0	0	0	0	0	0	3726	0	0	0	0	6599	0	0	0	0	6632	0	0	0	0	18735	0	0	0	57553
pos:67	100000	31.57	31.45	18.57	18.40	0.00	0	0	0	0	0	0	0	0	0	0	0	0	3328	0	0	0	0	0	0	0	0	0	2348	0	0	0	0	4546	0	0	0	0	4785	0	0	0	0	16716	0	0	0	68277
pos:68	100000	31.57	31.54	18.39	18.50	0.00	0	0	0	0	0	0	0	0	0	0	0	0	3051	0	0	0	0	0	0	0	0	0	1877	0	0	0	0	3364	0	0	0	0	5545	0	0	0	0	12746	0	0	0	73417
pos:69	100000	31.72	31.32	18.27	18.69	0.00	0	0	0	0	0	0	0	0	0	0	0	0	3155	0	0	0	0	0	0	0	0	0	1827	0	0	0	0	3281	0	0	0	0	5311	0	0	0	0	12095	0	0	0	74331
pos:70	100000	31.74	31.39	18.50	18.37	0.00	0	0	0	0	0	0	0	0	0	0	0	0	3564	0	0	0	0	0	0	0	0	0	1983	0	0	0	0	3408	0	0	0	0	5392	0	0	0	0	12306	0	0	0	73347
pos:71	100000	31.81	31.06	18.75	18.37	0.00	0	0	0	0	0	0	0	0	0	0	0	0	7178	0	0	0	0	0	0	0	0	0	3129	0	0	0	0	5093	0	0	0	0	7082	0	0	0	0	15502	0	0	0	62016
pos:72	100000	31.40	31.16	19.03	18.41	0.00	0	0	0	0	0	0	0	0	0	0	0	0	8939	0	0	0	0	0	0	0	0	0	3816	0	0	0	0	6159	0	0	0	0	8140	0	0	0	0	15056	0	0	0	57890
pos:73	100000	31.64	31.38	18.46	18.51	0.00	0	0	0	0	0	0	0	0	0	0	0	0	4036	0	0	0	0	0	0	0	0	0	2669	0	0	0	0	4927	0	0	0	0	7080	0	0	0	0	14006	0	0	0	67282
pos:74	100000	31.73	31.29	18.32	18.66	0.00	0	0	0	0	0	0	0	0	0	0	0	0	2072	0	0	0	0	0	0	0	0	0	1480	0	0	0	0	2838	0	0	0	0	4980	0	0	0	0	11838	0	0	0	76792
pos:75	100000	31.74	31.57	18.39	18.29	0.00	0	0	0	0	0	0	0	0	0	0	0	0	2666	0	0	0	0	0	0	0	0	0	1461	0	0	0	0	2608	0	0	0	0	4625	0	0	0	0	11448	0	0	0	77192
pos:76	100000	31.59	31.53	18.33	18.55	0.00	0	0	0	0	0	0	0	0	0	0	0	0	2163	0	0	0	0	0	0	0	0	0	2559	0	0	0	0	2581	0	0	0	0	6912	0	0	0	0	20066	0	0	0	65719
pos:77	100000	31.76	31.42	18.70	18.13	0.00	0	0	0	0	0	0	0	0	0	0	0	0	14255	0	0	0	0	0	0	0	0	0	4775	0	0	0	0	7677	0	0	0	0	9687	0	0	0	0	19117	0	0	0	44489
pos:78	100000	31.87	31.07	18.56	18.50	0.00	0	0	0	0	0	0	0	0	0	0	0	0	5762	0	0	0	0	0	0	0	0	0	5659	0	0	0	0	4096	0	0	0	0	9332	0	0	0	0	18377	0	0	0	56774
pos:79	100000	31.59	31.45	18.43	18.53	0.00	0	0	0	0	0	0	0	0	0	0	0	0	4372	0	0	0	0	0	0	0	0	0	4075	0	0	0	0	3155	0	0	0	0	7146	0	0	0	0	16129	0	0	0	65123
pos:80	100000	31.77	31.28	18.73	18.22	0.00	0	0	0	0	0	0	0	0	0	0	0	0	13269	0	0	0	0	0	0	0	0	0	6126	0	0	0	0	5747	0	0	0	0	8827	0	0	0	0	17421	0	0	0	48610
pos:81	100000	31.78	31.34	18.36	18.52	0.00	0	0	0	0	0	0	0	0	0	0	0	0	7813	0	0	0	0	0	0	0	0	0	6357	0	0	0	0	4583	0	0	0	0	9123	0	0	0	0	16825	0	0	0	55299
pos:82	100000	31.13	31.53	18.73	18.62	0.00	0	0	0	0	0	0	0	0	0	0	0	0	12438	0	0	0	0	0	0	0	0	0	6751	0	0	0	0	5301	0	0	0	0	9102	0	0	0	0	16603	0	0	0	49805
pos:83	100000	31.38	31.47	18.52	18.63	0.00	0	0	0	0	0	0	0	0	0	0	0	0	6126	0	0	0	0	0	0	0	0	0	6323	0	0	0	0	4289	0	0	0	0	8520	0	0	0	0	15931	0	0	0	58811
pos:84	100000	31.63	31.39	18.75	18.24	0.00	0	0	0	0	0	0	0	0	0	0	0	0	15973	0	0	0	0	0	0	0	0	0	7297	0	0	0	0	6028	0	0	0	0	9174	0	0	0	0	16464	0	0	0	45064
pos:85	100000	31.72	31.30	18.56	18.42	0.00	0	0	0	0	0	0	0	0	0	0	0	0	7158	0	0	0	0	0	0	0	0	0	7497	0	0	0	0	4814	0	0	0	0	10019	0	0	0	0	17466	0	0	0	53046
pos:86	100000	32.26	31.13	18.53	18.09	0.00	0	0	0	0	0	0	0	0	0	0	0	0	16460	0	0	0	0	0	0	0	0	0	7983	0	0	0	0	6354	0	0	0	0	9890	0	0	0	0	17557	0	0	0	41756
pos:87	100000	31.35	31.34	18.95	18.36	0.00	0	0	0	0	0	0	0	0	0	0	0	0	19635	0	0	0	0	0	0	0	0	0	9980	0	0	0	0	6645	0	0	0	0	10702	0	0	0	0	17116	0	0	0	35922
pos:88	100000	31.99	31.27	18.38	18.36	0.00	0	0	0	0	0	0	0	0	0	0	0	0	21238	0	0	0	0	0	0	0	0	0	11569	0	0	0	0	7140	0	0	0	0	11412	0	0	0	0	16996	0	0	0	31645
pos:89	100000	32.39	30.79	18.72	18.10	0.00	0	0	0	0	0	0	0	0	0	0	0	0	22077	0	0	0	0	0	0	0	0	0	12270	0	0	0	0	7411	0	0	0	0	11800	0	0	0	0	16978	0	0	0	29464
pos:90	100000	32.13	32.26	19.02	16.59	0.00	0	0	0	0	0	0	0	0	0	0	0	0	39004	0	0	0	0	0	0	0	0	0	13195	0	0	0	0	7687	0	0	0	0	9895	0	0	0	0	13042	0	0	0	17177
pos:91	100000	31.54	31.34	18.17	18.96	0.00	0	0	0	0	0	0	0	0	0	0	0	0	5981	0	0	0	0	0	0	0	0	0	12193	0	0	0	0	7805	0	0	0	0	15677	0	0	0	0	22443	0	0	0	35901
pos:92	100000	31.72	31.50	18.51	18.27	0.00	0	0	0	0	0	0	0	0	0	0	0	0	10435	0	0	0	0	0	0	0	0	0	8392	0	0	0	0	5894	0	0	0	0	12386	0	0	0	0	20632	0	0	0	42261
pos:93	100000	31.72	31.36	18.38	18.54	0.00	0	0	0	0	0	0	0	0	0	0	0	0	2930	0	0	0	0	0	0	0	0	0	4089	0	0	0	0	3097	0	0	0	0	8220	0	0	0	0	18246	0	0	0	63418
pos:94	100000	31.71	31.33	18.44	18.52	0.00	0	0	0	0	0	0	0	0	0	0	0	0	4536	0	0	0	0	0	0	0	0	0	3479	0	0	0	0	3189	0	0	0	0	6861	0	0	0	0	15980	0	0	0	65955
pos:95	100000	31.54	31.28	18.60	18.58	0.00	0	0	0	0	0	0	0	0	0	0	0	0	9994	0	0	0	0	0	0	0	0	0	5374	0	0	0	0	5223	0	0	0	0	8621	0	0	0	0	17514	0	0	0	53274
pos:96	100000	31.62	31.51	18.41	18.46	0.00	0	0	0	0	0	0	0	0	0	0	0	0	2765	0	0	0	0	0	0	0	0	0	3354	0	0	0	0	2647	0	0	0	0	6877	0	0	0	0	15408	0	0	0	68949
pos:97	100000	31.58	31.34	18.63	18.45	0.00	0	0	0	0	0	0	0	0	0	0	0	0	4457	0	0	0	0	0	0	0	0	0	3151	0	0	0	0	2945	0	0	0	0	6286	0	0	0	0	15211	0	0	0	67950
pos:98	100000	32.04	31.13	18.64	18.18	0.00	0	0	0	0	0	0	0	0	0	0	0	0	14703	0	0	0	0	0	0	0	0	0	6202	0	0	0	0	6015	0	0	0	0	8723	0	0	0	0	17402	0	0	0	46955
pos:99	100000	31.61	31.05	18.84	18.51	0.00	0	0	0	0	0	0	0	0	0	0	0	0	10396	0	0	0	0	0	0	0	0	0	7483	0	0	0	0	5294	0	0	0	0	10167	0	0	0	0	17949	0	0	0	48711
pos:100	100000	31.86	31.15	18.41	18.58	0.00	0	0	0	0	0	0	0	0	0	0	0	0	10276	0	0	0	0	0	0	0	0	0	7207	0	0	0	0	5077	0	0	0	0	10777	0	0	0	0	16155	0	0	0	50508
pos:101	100000	31.45	31.48	18.75	18.33	0.00	0	0	0	0	0	0	0	0	0	0	0	0	20555	0	0	0	0	0	0	0	0	0	9192	0	0	0	0	7096	0	0	0	0	10310	0	0	0	0	17405	0	0	0	35442
pos:102	100000	31.21	31.21	18.71	18.86	0.00	0	0	0	0	0	0	0	0	0	0	0	0	12572	0	0	0	0	0	0	0	0	0	10050	0	0	0	0	6107	0	0	0	0	12788	0	0	0	0	17152	0	0	0	41331
pos:103	100000	32.43	31.80	18.17	17.60	0.00	0	0	0	0	0	0	0	0	0	0	0	0	32153	0	0	0	0	0	0	0	0	0	11055	0	0	0	0	7838	0	0	0	0	10451	0	0	0	0	15427	0	0	0	23076
pos:104	100000	31.44	31.49	18.67	18.39	0.00	0	0	0	0	0	0	0	0	0	0	0	0	10041	0	0	0	0	0	0	0	0	0	12066	0	0	0	0	6739	0	0	0	0	15561	0	0	0	0	18870	0	0	0	36723
pos:105	100000	31.16	31.18	18.74	18.91	0.00	0	0	0	0	0	0	0	0	0	0	0	0	5668	0	0	0	0	0	0	0	0	0	6815	0	0	0	0	4750	0	0	0	0	12046	0	0	0	0	19296	0	0	0	51425
pos:106	100000	31.54	31.29	18.83	18.34	0.00	0	0	0	0	0	0	0	0	0	0	0	0	10947	0	0	0	0	0	0	0	0	0	6561	0	0	0	0	5322	0	0	0	0	10893	0	0	0	0	16990	0	0	0	49287
pos:107	100000	32.02	31.09	18.26	18.63	0.00	0	0	0	0	0	0	0	0	0	0	0	0	7323	0	0	0	0	0	0	0	0	0	6321	0	0	0	0	4991	0	0	0	0	10551	0	0	0	0	16599	0	0	0	54215
pos:108	100000	31.97	31.20	18.71	18.12	0.00	0	0	0	0	0	0	0	0	0	0	0	0	15467	0	0	0	0	0	0	0	0	0	7478	0	0	0	0	6804	0	0	0	0	10951	0	0	0	0	16652	0	0	0	42648
pos:109	100000	31.45	31.23	18.70	18.61	0.00	0	0	0	0	0	0	0	0	0	0	0	0	10368	0	0	0	0	0	0	0	0	0	8098	0	0	0	0	5756	0	0	0	0	11374	0	0	0	0	16344	0	0	0	48060
pos:110	100000	31.44	31.59	18.36	18.62	0.00	0	0	0	0	0	0	0	0	0	0	0	0	6873	0	0	0	0	0	0	0	0	0	6568	0	0	0	0	4787	0	0	0	0	10330	0	0	0	0	16039	0	0	0	55403
pos:111	100000	31.44	31.27	18.71	18.58	0.00	0	0	0	0	0	0	0	0	0	0	0	0	8366	0	0	0	0	0	0	0	0	0	6061	0	0	0	0	5089	0	0	0	0	9706	0	0	0	0	16230	0	0	0	54548
pos:112	100000	31.28	31.73	18.93	18.06	0.00	0	0	0	0	0	0	0	0	0	0	0	0	21795	0	0	0	0	0	0	0	0	0	8057	0	0	0	0	6608	0	0	0	0	9637	0	0	0	0	14651	0	0	0	39252
pos:113	100000	31.66	31.38	19.16	17.80	0.00	0	0	0	0	0	0	0	0	0	0	0	0	25734	0	0	0	0	0	0	0	0	0	11588	0	0	0	0	7560	0	0	0	0	11143	0	0	0	0	14025	0	0	0	29950
pos:114	100000	31.71	31.48	18.42	18.39	0.00	0	0	0	0	0	0	0	0	0	0	0	0	10451	0	0	0	0	0	0	0	0	0	15084	0	0	0	0	3738	0	0	0	0	13527	0	0	0	0	20446	0	0	0	36754
pos:115	100000	31.72	31.24	18.47	18.57	0.00	0	0	0	0	0	0	0	0	0	0	0	0	3186	0	0	0	0	0	0	0	0	0	6807	0	0	0	0	2374	0	0	0	0	10979	0	0	0	0	24678	0	0	0	51976
pos:116	100000	31.52	31.57	18.93	17.98	0.00	0	0	0	0	0	0	0	0	0	0	0	0	24462	0	0	0	0	0	0	0	0	0	9217	0	0	0	0	7246	0	0	0	0	11288	0	0	0	0	18801	0	0	0	28986
pos:117	100000	31.30	31.43	19.30	17.96	0.00	0	0	0	0	0	0	0	0	0	0	0	0	28409	0	0	0	0	0	0	0	0	0	13515	0	0	0	0	5815	0	0	0	0	11963	0	0	0	0	16773	0	0	0	23525
pos:118	100000	31.46	30.86	18.64	19.03	0.00	0	0	0	0	0	0	0	0	0	0	0	0	18936	0	0	0	0	0	0	0	0	0	16310	0	0	0	0	4522	0	0	0	0	12883	0	0	0	0	18662	0	0	0	28687
pos:119	100000	31.40	31.09	18.91	18.60	0.00	0	0	0	0	0	0	0	0	0	0	0	0	16539	0	0	0	0	0	0	0	0	0	15047	0	0	0	0	4674	0	0	0	0	12694	0	0	0	0	19605	0	0	0	31441
pos:120	100000	31.62	31.08	18.63	18.67	0.00	0	0	0	0	0	0	0	0	0	0	0	0	15017	0	0	0	0	0	0	0	0	0	13437	0	0	0	0	4644	0	0	0	0	12275	0	0	0	0	19805	0	0	0	34822
pos:121	100000	31.68	31.05	18.54	18.73	0.00	0	0	0	0	0	0	0	0	0	0	0	0	11180	0	0	0	0	0	0	0	0	0	11886	0	0	0	0	5929	0	0	0	0	9964	0	0	0	0	21844	0	0	0	39197
pos:122	100000	31.81	31.00	18.56	18.63	0.00	0	0	0	0	0	0	0	0	0	0	0	0	10737	0	0	0	0	0	0	0	0	0	9933	0	0	0	0	6047	0	0	0	0	9663	0	0	0	0	21990	0	0	0	41630
pos:123	100000	31.90	31.17	18.53	18.40	0.00	0	0	0	0	0	0	0	0	0	0	0	0	13478	0	0	0	0	0	0	0	0	0	9670	0	0	0	0	6990	0	0	0	0	9658	0	0	0	0	21558	0	0	0	38646
pos:124	100000	31.62	31.22	18.58	18.58	0.00	0	0	0	0	0	0	0	0	0	0	0	0	6838	0	0	0	0	0	0	0	0	0	8732	0	0	0	0	5042	0	0	0	0	9548	0	0	0	0	23134	0	0	0	46706
pos:125	100000	31.81	31.06	18.61	18.51	0.00	0	0	0	0	0	0	0	0	0	0	0	0	8964	0	0	0	0	0	0	0	0	0	5352	0	0	0	0	5016	0	0	0	0	8659	0	0	0	0	23651	0	0	0	48358
pos:126	100000	31.59	31.39	18.56	18.46	0.00	0	0	0	0	0	0	0	0	0	0	0	0	12997	0	0	0	0	0	0	0	0	0	5858	0	0	0	0	6012	0	0	0	0	9006	0	0	0	0	22921	0	0	0	43206
pos:127	100000	31.53	31.53	18.55	18.39	0.00	0	0	0	0	0	0	0	0	0	0	0	0	15367	0	0	0	0	0	0	0	0	0	6639	0	0	0	0	5906	0	0	0	0	9110	0	0	0	0	21858	0	0	0	41120
pos:128	100000	31.77	31.13	18.63	18.46	0.00	0	0	0	0	0	0	0	0	0	0	0	0	16501	0	0	0	0	0	0	0	0	0	7601	0	0	0	0	6424	0	0	0	0	9856	0	0	0	0	22415	0	0	0	37203
pos:129	100000	31.47	31.11	18.59	18.83	0.00	0	0	0	0	0	0	0	0	0	0	0	0	13863	0	0	0	0	0	0	0	0	0	7507	0	0	0	0	5965	0	0	0	0	9728	0	0	0	0	22879	0	0	0	40058
pos:130	100000	31.59	31.15	18.40	18.86	0.00	0	0	0	0	0	0	0	0	2594	0	0	0	15205	0	0	0	0	0	0	0	0	0	8217	0	0	0	0	5987	0	0	0	0	13758	0	0	0	0	19098	0	0	0	35141
pos:131	100000	31.42	31.75	18.91	17.92	0.00	0	0	0	0	0	0	0	0	0	0	0	0	28683	0	0	0	0	0	0	0	0	0	8470	0	0	0	0	7534	0	0	0	0	9562	0	0	0	0	19782	0	0	0	25969
pos:132	100000	31.55	31.26	18.72	18.46	0.00	0	0	0	0	0	0	0	0	0	0	0	0	12330	0	0	0	0	0	0	0	0	0	11125	0	0	0	0	5616	0	0	0	0	11489	0	0	0	0	23930	0	0	0	35510
pos:133	100000	31.51	31.14	18.63	18.71	0.00	0	0	0	0	0	0	0	0	984	0	0	0	5282	0	0	0	0	0	0	0	0	0	6490	0	0	0	0	4219	0	0	0	0	12863	0	0	0	0	25076	0	0	0	45086
pos:134	100000	31.72	31.21	18.49	18.57	0.00	0	0	0	0	0	0	0	0	1560	0	0	0	8631	0	0	0	0	0	0	0	0	0	5648	0	0	0	0	5042	0	0	0	0	13018	0	0	0	0	21217	0	0	0	44884
pos:135	100000	31.50	31.43	18.65	18.41	0.00	0	0	0	0	0	0	0	0	1944	0	0	0	11485	0	0	0	0	0	0	0	0	0	6832	0	0	0	0	5579	0	0	0	0	14173	0	0	0	0	20350	0	0	0	39637
pos:136	100000	31.68	31.69	18.91	17.72	0.00	0	0	0	0	0	0	0	0	0	0	0	0	30597	0	0	0	0	0	0	0	0	0	8141	0	0	0	0	8311	0	0	0	0	9940	0	0	0	0	20183	0	0	0	22828
pos:137	100000	31.75	31.44	18.37	18.44	0.00	0	0	0	0	0	0	0	0	2089	0	0	0	13305	0	0	0	0	0	0	0	0	0	11230	0	0	0	0	6034	0	0	0	0	14886	0	0	0	0	22160	0	0	0	30296
pos:138	100000	31.84	31.38	18.44	18.33	0.00	0	0	0	0	0	0	0	0	1712	0	0	0	10273	0	0	0	0	0	0	0	0	0	8803	0	0	0	0	5606	0	0	0	0	14069	0	0	0	0	22862	0	0	0	36675
pos:139	100000	31.20	31.50	18.62	18.67	0.00	0	0	0	0	0	0	0	0	1208	0	0	0	7651	0	0	0	0	0	0	0	0	0	6629	0	0	0	0	4776	0	0	0	0	13027	0	0	0	0	22664	0	0	0	44045
pos:140	100000	31.60	31.23	18.48	18.70	0.00	0	0	0	0	0	0	0	0	1422	0	0	0	8094	0	0	0	0	0	0	0	0	0	5684	0	0	0	0	4847	0	0	0	0	13174	0	0	0	0	21664	0	0	0	45115
pos:141	100000	31.41	31.51	18.64	18.45	0.00	0	0	0	0	0	0	0	0	1477	0	0	0	8818	0	0	0	0	0	0	0	0	0	6014	0	0	0	0	4974	0	0	0	0	13359	0	0	0	0	21123	0	0	0	44235
pos:142	100000	31.49	31.56	18.75	18.21	0.00	0	0	0	0	0	0	0	0	2952	0	0	0	16987	0	0	0	0	0	0	0	0	0	8302	0	0	0	0	6594	0	0	0	0	14891	0	0	0	0	18696	0	0	0	31578
pos:143	100000	31.46	31.23	18.62	18.68	0.00	0	0	0	0	0	0	0	0	1763	0	0	0	10482	0	0	0	0	0	0	0	0	0	8216	0	0	0	0	5364	0	0	0	0	14222	0	0	0	0	22659	0	0	0	37294
pos:144	100000	31.05	31.43	18.86	18.65	0.00	0	0	0	0	0	0	0	0	4104	0	0	0	23695	0	0	0	0	0	0	0	0	0	9489	0	0	0	0	6837	0	0	0	0	14349	0	0	0	0	17157	0	0	0	24369
pos:145	100000	31.44	31.23	18.44	18.89	0.00	0	0	0	0	0	0	0	0	1677	0	0	0	11376	0	0	0	0	0	0	0	0	0	10428	0	0	0	0	10101	0	0	0	0	16561	0	0	0	0	16239	0	0	0	33618
pos:146	100000	31.41	31.08	19.04	18.47	0.00	0	0	0	0	0	0	0	0	3123	0	0	0	19448	0	0	0	0	0	0	0	0	0	9116	0	0	0	0	9602	0	0	0	0	14651	0	0	0	0	14824	0	0	0	29236
pos:147	100000	31.11	31.37	19.04	18.48	0.00	0	0	0	0	0	0	0	0	4744	0	0	0	27333	0	0	0	0	0	0	0	0	0	10959	0	0	0	0	7053	0	0	0	0	14308	0	0	0	0	16027	0	0	0	19576
pos:148	100000	31.38	31.05	19.17	18.40	0.00	0	0	0	0	0	0	0	0	4487	0	0	0	26753	0	0	0	0	0	0	0	0	0	13559	0	0	0	0	6937	0	0	0	0	14102	0	0	0	0	16424	0	0	0	17738
pos:149	100000	31.30	31.05	18.83	18.82	0.00	0	0	0	0	0	0	0	0	2313	0	0	0	16224	0	0	0	0	0	0	0	0	0	13607	0	0	0	0	10976	0	0	0	0	16760	0	0	0	0	14934	0	0	0	25186
pos:150	100000	31.38	31.04	18.79	18.79	0.00	0	0	0	0	0	0	0	0	2439	0	0	0	15943	0	0	0	0	0	0	0	0	0	10277	0	0	0	0	10001	0	0	0	0	16133	0	0	0	0	15672	0	0	0	29535
```
output file: `summary.log`
```
Total read number:      100000
Total base number:      15000000
Max read length:        150
Reads average length:   150.0
Total GC content:       36.84
```

## Task list

*   [x] add base distrbution plot

*   [x] add output matrix file and summary log

*   [x] add base quality plot

*   [ ] add command line args

*   [ ] rename output png file

** any bugs please report issues **💖 
