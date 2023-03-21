#!/usr/bin/env perl
#pk

##########################################################################################
#	  This file is part of proteinortho.
#	  (C) 2009 Marcus Lechner
# 
#	  proteinortho is free software; you can redistribute it and/or modify
#	  it under the terms of the GNU General Public License as published
#	  by the Free Software Foundation; either version 2, or (at your
#	  option) any later version.
#
#	  proteinortho is distributed in the hope that it will be useful, but
#	  WITHOUT ANY WARRANTY; without even the implied warranty of
#	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#	  General Public License for more details.
#
#	  You should have received a copy of the GNU General Public License
#	  along with proteinortho; see the file COPYING.  If not, write to the
#	  Free Software Foundation, Inc., 59 Temple Place - Suite 330,
#	  Boston, MA 02111-1307, USA.	
##########################################################################################

##########################################################################################
# About
##########################################################################################
# proteinortho2html.pl
# input .proteinortho.tsv file and (optional) the used genomes
# output HTML formatted proteinortho.tsv
#
# @author Paul Klemm
# @email klemmp@staff.uni-marburg.de
# @company Bioinformatics, University of Leipzig
# @version 2
# @date 8-12-2020
#
##########################################################################################

if(scalar(@ARGV)<1){
	print STDERR "USAGE: proteinortho2html.pl <myproject.proteinortho> (<fasta1> <fasta2> ...)\nthe first argument points to the proteinortho output (tsv)-file. Any further (optional) files should be fasta files, for conversion of the identifier to a proper gene name/ describtion. The HTML output is printed to stdout, use '>' to write the html output to a file.\n";
	exit 1;
}elsif($ARGV[0] =~m/-he?l?p?$/){
	print STDERR "USAGE: proteinortho2html.pl <myproject.proteinortho> (<fasta1> <fasta2> ...)\nthe first argument points to the proteinortho output (tsv)-file. Any further (optional) files should be fasta files, for conversion of the identifier to a proper gene name/ describtion. The HTML output is printed to stdout, use '>' to write the html output to a file.\n";
	exit 0;
}

#
## get the id conversion data (id->gene describtion)
#

my %mappingdata;
my $proteinorthotsvfile;
	
for(my $argi=0;$argi<scalar(@ARGV);$argi++){
	if($ARGV[$argi]=~m/\.tsv/){$proteinorthotsvfile=$ARGV[$argi];next;}
	if($ARGV[$argi]=~m/proteinortho2html/){next;}
	open(FILE,"<",$ARGV[$argi]);
	while(<FILE>){
		if(substr($_,0,1) ne ">"){next;}
		chomp;
		my $curline=$_;
		if($curline=~m/^>([^ ]+) ([^=]+) [^ ]+=/){
			$mappingdata{$1}=$2;
			$mappingdata{$1}=~s/^([A-Z][a-z]+ [a-z]+)//;
			$mappingdata{$1}=~s/^([^ ]+ )([A-Z][a-z]+ [a-z]+)/$1/;
		}elsif($curline=~m/^>([^ ]+) .*GN=([^ ]+)/){
			$mappingdata{$1}=$2;
			$mappingdata{$1}=~s/^([A-Z][a-z]+ [a-z]+)//;
			$mappingdata{$1}=~s/^([^ ]+ )([A-Z][a-z]+ [a-z]+)/$1/;
		}elsif($curline=~m/^>([^ ]+) ([^=()]+)/){
			$mappingdata{$1}=$2;
			$mappingdata{$1}=~s/^([A-Z][a-z]+ [a-z]+)//;
			$mappingdata{$1}=~s/^([^ ]+ )([A-Z][a-z]+ [a-z]+)/$1/;
		}
	}
	close(FILE);
}

if(!-e $proteinorthotsvfile){
	print STDERR "USAGE: proteinortho2html.pl <myproject.proteinortho> (<fasta1> <fasta2> ...)\nthe first argument points to the proteinortho output (tsv)-file. Any further (optional) files should be fasta files, for conversion of the identifier to a proper gene name/ describtion. The HTML output is printed to stdout, use '>' to write the html output to a file.\n";
	exit 1;
}

#
## HEAD STUFF
#


print '<!DOCTYPE HTML>
<html>
  <head>
    <meta charset="UTF-8">
    <title>Proteinortho 6</title>
	<script type="text/javascript">

		// ***** clusterize.min.js ***** (source:https://github.com/NeXTs/Clusterize.js)
;(function(q,n){"undefined"!=typeof module?module.exports=n():"function"==typeof define&&"object"==typeof define.amd?define(n):this[q]=n()})("Clusterize",function(){function q(b,a,c){return a.addEventListener?a.addEventListener(b,c,!1):a.attachEvent("on"+b,c)}function n(b,a,c){return a.removeEventListener?a.removeEventListener(b,c,!1):a.detachEvent("on"+b,c)}function r(b){return"[object Array]"===Object.prototype.toString.call(b)}function m(b,a){return window.getComputedStyle?window.getComputedStyle(a)[b]:
a.currentStyle[b]}var l=function(){for(var b=3,a=document.createElement("b"),c=a.all||[];a.innerHTML="\x3c!--[if gt IE "+ ++b+"]><i><![endif]--\x3e",c[0];);return 4<b?b:document.documentMode}(),x=navigator.platform.toLowerCase().indexOf("mac")+1,p=function(b){if(!(this instanceof p))return new p(b);var a=this,c={rows_in_block:50,blocks_in_cluster:4,tag:null,show_no_data_row:!0,no_data_class:"clusterize-no-data",no_data_text:"No data",keep_parity:!0,callbacks:{}};a.options={};for(var d="rows_in_block blocks_in_cluster show_no_data_row no_data_class no_data_text keep_parity tag callbacks".split(" "),
f=0,h;h=d[f];f++)a.options[h]="undefined"!=typeof b[h]&&null!=b[h]?b[h]:c[h];c=["scroll","content"];for(f=0;d=c[f];f++)if(a[d+"_elem"]=b[d+"Id"]?document.getElementById(b[d+"Id"]):b[d+"Elem"],!a[d+"_elem"])throw Error("Error! Could not find "+d+" element");a.content_elem.hasAttribute("tabindex")||a.content_elem.setAttribute("tabindex",0);var e=r(b.rows)?b.rows:a.fetchMarkup(),g={};b=a.scroll_elem.scrollTop;a.insertToDOM(e,g);a.scroll_elem.scrollTop=b;var k=!1,m=0,l=!1,t=function(){x&&(l||(a.content_elem.style.pointerEvents=
"none"),l=!0,clearTimeout(m),m=setTimeout(function(){a.content_elem.style.pointerEvents="auto";l=!1},50));k!=(k=a.getClusterNum())&&a.insertToDOM(e,g);a.options.callbacks.scrollingProgress&&a.options.callbacks.scrollingProgress(a.getScrollProgress())},u=0,v=function(){clearTimeout(u);u=setTimeout(a.refresh,100)};q("scroll",a.scroll_elem,t);q("resize",window,v);a.destroy=function(b){n("scroll",a.scroll_elem,t);n("resize",window,v);a.html((b?a.generateEmptyRow():e).join(""))};a.refresh=function(b){(a.getRowsHeight(e)||
b)&&a.update(e)};a.update=function(b){e=r(b)?b:[];b=a.scroll_elem.scrollTop;e.length*a.options.item_height<b&&(k=a.scroll_elem.scrollTop=0);a.insertToDOM(e,g);a.scroll_elem.scrollTop=b};a.clear=function(){a.update([])};a.getRowsAmount=function(){return e.length};a.getScrollProgress=function(){return this.options.scroll_top/(e.length*this.options.item_height)*100||0};var w=function(b,c){var d=r(c)?c:[];d.length&&(e="append"==b?e.concat(d):d.concat(e),a.insertToDOM(e,g))};a.append=function(a){w("append",
a)};a.prepend=function(a){w("prepend",a)}};p.prototype={constructor:p,fetchMarkup:function(){for(var b=[],a=this.getChildNodes(this.content_elem);a.length;)b.push(a.shift().outerHTML);return b},exploreEnvironment:function(b,a){var c=this.options;c.content_tag=this.content_elem.tagName.toLowerCase();b.length&&(l&&9>=l&&!c.tag&&(c.tag=b[0].match(/<([^>\s/]*)/)[1].toLowerCase()),1>=this.content_elem.children.length&&(a.data=this.html(b[0]+b[0]+b[0])),c.tag||(c.tag=this.content_elem.children[0].tagName.toLowerCase()),
this.getRowsHeight(b))},getRowsHeight:function(b){var a=this.options,c=a.item_height;a.cluster_height=0;if(b.length&&(b=this.content_elem.children,b.length)){var d=b[Math.floor(b.length/2)];a.item_height=d.offsetHeight;"tr"==a.tag&&"collapse"!=m("borderCollapse",this.content_elem)&&(a.item_height+=parseInt(m("borderSpacing",this.content_elem),10)||0);"tr"!=a.tag&&(b=parseInt(m("marginTop",d),10)||0,d=parseInt(m("marginBottom",d),10)||0,a.item_height+=Math.max(b,d));a.block_height=a.item_height*a.rows_in_block;
a.rows_in_cluster=a.blocks_in_cluster*a.rows_in_block;a.cluster_height=a.blocks_in_cluster*a.block_height;return c!=a.item_height}},getClusterNum:function(){this.options.scroll_top=this.scroll_elem.scrollTop;return Math.floor(this.options.scroll_top/(this.options.cluster_height-this.options.block_height))||0},generateEmptyRow:function(){var b=this.options;if(!b.tag||!b.show_no_data_row)return[];var a=document.createElement(b.tag),c=document.createTextNode(b.no_data_text);a.className=b.no_data_class;
if("tr"==b.tag){var d=document.createElement("td");d.colSpan=100;d.appendChild(c)}a.appendChild(d||c);return[a.outerHTML]},generate:function(b,a){var c=this.options,d=b.length;if(d<c.rows_in_block)return{top_offset:0,bottom_offset:0,rows_above:0,rows:d?b:this.generateEmptyRow()};var f=Math.max((c.rows_in_cluster-c.rows_in_block)*a,0),h=f+c.rows_in_cluster,e=Math.max(f*c.item_height,0);c=Math.max((d-h)*c.item_height,0);d=[];var g=f;for(1>e&&g++;f<h;f++)b[f]&&d.push(b[f]);return{top_offset:e,bottom_offset:c,
rows_above:g,rows:d}},renderExtraTag:function(b,a){var c=document.createElement(this.options.tag);c.className=["clusterize-extra-row","clusterize-"+b].join(" ");a&&(c.style.height=a+"px");return c.outerHTML},insertToDOM:function(b,a){this.options.cluster_height||this.exploreEnvironment(b,a);var c=this.generate(b,this.getClusterNum()),d=c.rows.join(""),f=this.checkChanges("data",d,a),h=this.checkChanges("top",c.top_offset,a),e=this.checkChanges("bottom",c.bottom_offset,a),g=this.options.callbacks,
k=[];f||h?(c.top_offset&&(this.options.keep_parity&&k.push(this.renderExtraTag("keep-parity")),k.push(this.renderExtraTag("top-space",c.top_offset))),k.push(d),c.bottom_offset&&k.push(this.renderExtraTag("bottom-space",c.bottom_offset)),g.clusterWillChange&&g.clusterWillChange(),this.html(k.join("")),"ol"==this.options.content_tag&&this.content_elem.setAttribute("start",c.rows_above),this.content_elem.style["counter-increment"]="clusterize-counter "+(c.rows_above-1),g.clusterChanged&&g.clusterChanged()):
e&&(this.content_elem.lastChild.style.height=c.bottom_offset+"px")},html:function(b){var a=this.content_elem;if(l&&9>=l&&"tr"==this.options.tag){var c=document.createElement("div");for(c.innerHTML="<table><tbody>"+b+"</tbody></table>";b=a.lastChild;)a.removeChild(b);for(c=this.getChildNodes(c.firstChild.firstChild);c.length;)a.appendChild(c.shift())}else a.innerHTML=b},getChildNodes:function(b){b=b.children;for(var a=[],c=0,d=b.length;c<d;c++)a.push(b[c]);return a},checkChanges:function(b,a,c){var d=
a!=c[b];c[b]=a;return d}};return p});
		// ***** clusterize.min.js *****

		var clusterize;
		
		function isEmpty() {
			input = document.getElementById("searchField");
			if(input.value==""){
				document.getElementById("countspanA").textContent="nothing";
				document.getElementById("countspanB").textContent=data_raw.length+"/"+data_raw.length;
				data_filtered_idx=[];
				for (var i = 0; i < data_raw.length; i++) {
					data_filtered_idx.push(i);
				}
				updateTable();
			}
		}
		function extractContent(s) {
			var span = document.createElement(\'span\');
			span.innerHTML = s;

			txtValue = span.textContent.toUpperCase() || span.innerText.toUpperCase();
			txtValue = txtValue.replace(\'(+)\',\'\');
			txtValue = txtValue.replace(\'(-)\',\'\');
			txtValue = txtValue.replace(\' \',\'\');
			return span.textContent.toUpperCase() || span.innerText.toUpperCase();
		};

		var data_filtered_idx=[];

		function filterTable(typeofsearch=0) {

			updateTable();
			
			data_filtered_idx=[];
			input = document.getElementById("searchField");
			filter = input.value.toUpperCase();

			if(useMappedData){
				for (var i = 0; i < data_raw_mapped.length; i++) {
					tds = data_raw_mapped[i].split(/<\/td><td[^>]*>/);
					
					for (j = 3; j < tds.length; j++) {
						if(selectedCols.indexOf(j)==-1){continue;}

						txtValue=extractContent(tds[j]);

						if(typeofsearch==0){
							if (txtValue.indexOf(filter) != -1) {
								data_filtered_idx.push(i);
								break;
							}
						}else{
							var words = txtValue.split(/[,]/);
							if (words.includes(filter)) {
								data_filtered_idx.push(i);
								break;
							}
						}
					}
				}
			}else{
				for (var i = 0; i < data_raw.length; i++) {
					tds = data_raw[i].split(/<\/td><td[^>]*>/);
					
					for (j = 3; j < tds.length; j++) {
						if(selectedCols.indexOf(j)==-1){continue;}

						txtValue=extractContent(tds[j]);

						if(typeofsearch==0){
							if (txtValue.indexOf(filter) != -1) {
								data_filtered_idx.push(i);
								break;
							}
						}else{
							var words = txtValue.split(/[,]/);
							if (words.includes(filter)) {
								data_filtered_idx.push(i);
								break;
							}
						}
					}
				}
			}
			
			if(typeofsearch==0){
				document.getElementById("countspanA").textContent="*"+filter+"*";
				document.getElementById("countspanB").textContent=data_filtered_idx.length+"/"+data_raw.length;
			}else{
				document.getElementById("countspanA").textContent="\'"+filter+"\'";
				document.getElementById("countspanB").textContent=data_filtered_idx.length+"/"+data_raw.length;
			}

			updateTable();
		}

		function showpopup(id) {
			var popup = document.getElementById("myPopupText"+id);
			if(popup.style.display=="inline-block"){
				popup.style.display="none";
			}else{
				popup.style.display="inline-block";
			}

			var popupbutton = document.getElementById("myPopup"+id);
			if(popupbutton.innerHTML.includes(\'(+)\')){
				popupbutton.innerHTML = popupbutton.innerHTML.replace( /\(\+\)/,\'(-)\');
			}else{
				 popupbutton.innerHTML = popupbutton.innerHTML.replace( /\(\-\)/,\'(+)\');
			}
		}
	</script>
    <style>
		table.greyGridTable {
		  border: 2px solid #FFFFFF;
/*		  width: 100%; */
		  text-align: center;
		  min-width:100%;
		  border-collapse: collapse;
		  table-layout: fixed;
		  z-index: 999;
		}
		table.greyGridTable td, table.greyGridTable th {
		  border: 1px solid #FFFFFF;
		  padding: 3px 4px;
		}
		table.greyGridTable tbody td {
		  font-size: 13px;
		  word-wrap: break-word;         /* All browsers since IE 5.5+ */
		    overflow-wrap: break-word;     /* Renamed property in CSS3 draft spec */
		}
		table.greyGridTable thead {
		  background: #EBEBEB;

		  font-size: 15px;
		  font-weight: bold;
		  color: #333333;
		  text-align: center;

		    position: -webkit-sticky;
		    position: -moz-sticky; 
		    position: -ms-sticky;
		    position: -o-sticky;
		    position: sticky;
    		top: 0;
    		z-index: 99;

		}
		
		table.greyGridTable thead th:not(:nth-child(1)):not(:nth-child(2)):not(:nth-child(3)) {
	      
	      -webkit-transform: rotate(-0deg) translatex(-0%);
	      -moz-transform: rotate(-0deg) translatex(-0%);
	      -ms-transform: rotate(-0deg) translatex(-0%);
	      -o-transform: rotate(-0deg) translatex(-0%);
	      transform: rotate(-0deg) translatex(-0%);
	      font-weight: normal;
/*		  word-wrap: break-word;         
		    overflow-wrap: break-word;      */

		}
		table.greyGridTable thead th:nth-child(-n+2) {
	      
		  font-size:8px;
		  white-space: nowrap;

		}

		table.greyGridTable thead th:nth-child(-n+2) (height > 10px) {

	      -webkit-transform: rotate(-90deg) translatex(-0%);
	      -moz-transform: rotate(-90deg) translatex(-0%);
	      -ms-transform: rotate(-90deg) translatex(-0%);
	      -o-transform: rotate(-90deg) translatex(-0%);
	      transform: rotate(-90deg) translatex(-0%);

		}

		table.greyGridTable thead th:nth-child(-n+3) {
	      width:50px;
		  font-size:8px;
		}
		table.greyGridTable tfoot th:nth-child(-n+3) {
	      width:50px;
		  font-size:9px;
		}
		table.greyGridTable tbody th:nth-child(-n+3) {
	      width:50px;
		  font-size:9px;
		}

		table.greyGridTable tr td:nth-child(-n+3) {
			      
		  font-size:12px;
		  white-space: nowrap;

		}

		table.greyGridTable tfoot th {
		  position: -webkit-sticky;
		  position: sticky;
		    bottom: 0;
		  background: #EBEBEB;
		}

		table.greyGridTable tfoot {
		  font-size: 14px;
		  font-weight: bold;
		  color: #333333;
		}
		table.greyGridTable tfoot td {
		  font-size: 14px;
		}
		table.greyGridTable tr:hover td:not(:nth-child(-n+3)){
		    background: #ddd;
		}
		table.greyGridTable tr:hover td:nth-child(-n+3){
		    background: rgb(140, 184, 255,0.65);
		    cursor:pointer;
		}
		h3{
		  white-space: nowrap;
		  height:25px;
		}
		.popup {
		  position: relative;
		  display: inline-block;
		  cursor: pointer;
		  -webkit-user-select: none;
		  -moz-user-select: none;
		  -ms-user-select: none;
		  user-select: none;
		}
		.popuptext {
		  display: none;
		}
		.clusterize-scroll{
		    position: absolute;
		    top: 155px;
		    width: 100%;
		    left: 0;
		    right: 0;
		    padding:0;
		    margin:0;
		    bottom: 0;
		    max-height: none; /* <-- important one  */
		  	overflow: auto;
		}
		.clusterize-extra-row{
		  margin-top: 0 !important;
		  margin-bottom: 0 !important;
		}
		.clusterize-extra-row.clusterize-keep-parity{
		  display: none;
		}

		.clusterize-content{
		  outline: 0;
		  counter-reset: clusterize-counter;
		}
		.clusterize-no-data td{
		  text-align: center;
		}
		.overlay{
		  padding:5px;
		  background: #f0f2ff;
		  position:fixed;
		  width:100%;
		  height:100px;
		  left:0px;
		  border-top: 3px #dbdbdb solid;
		  border-bottom: 3px #dbdbdb solid;
		}
		.left{
			width: 40%;
		    background: rgb(255,255,255,0);
		    float: left;	
		    height:100%;
		}
		.middleA {    
			width: 10%;
		    height:90%;
		    float: left;	
		}
		.middleB {
		    float: left;	
		    height:90%;
		    margin-left:5%;
			width: 45%;
		}
		.middleA input{
			width:100%
		}
		.middleA input[type="checkbox"]{
			width:10px;
		}
		.sidewaysdiv {
			border-radius: 5px 5px 5px 5px;
	        user-select: none;
	        -moz-user-select: none;
	        -khtml-user-select: none;
	        -webkit-user-select: none;
	        -o-user-select: none;

		    /* Rotate from top left corner (not default) */
		    -webkit-transform-origin: 0 0;
		    -moz-transform-origin:    0 0;
		    -ms-transform-origin:     0 0;
		    -o-transform-origin:      0 0;
		    
		    -webkit-transform: rotate(90deg); 
		    -moz-transform:    rotate(90deg); 
		    -ms-transform:     rotate(90deg); 
		    -o-transform:      rotate(90deg); 
			  top:0px;
			  left:-5px;
			  font-size:14px;
			  position: relative;
			  
		    background: white;
		    color: black;
		    padding: 4px 5px;
		    margin: 0 0 3px 0;
		    width:80px;
		    height:16px;

		  	text-align: center;
		}
		.right {
			width: 150px;
			float: right;
		    background: rgb(255,255,255,0);
		}
	    select {
	      border: none;  box-shadow: 0px 2px 5px 1px rgba(0,0,0,0.3);
		  border-radius: 3px;
		  border-top-left-radius:0px;
		  border-bottom-left-radius:0px;
		  top:5px;
		  position: absolute;
	    }


	    input {
	      border: none;  box-shadow: 2px 2px 5px 1px rgba(0,0,0,0.3);
	  	  border-radius: 3px;
	    } 
	    h3{
	    	margin:4px;
	    	z-index:9999;
	    }
	    span{
	    	font-size:13px;
	    }
	    .modal {
		  display: none; /* Hidden by default */
		  position: fixed; /* Stay in place */
		  z-index: 1; /* Sit on top */
		  left: 0;
		  top: 0;
		  width: 100%; /* Full width */
		  height: 100%; /* Full height */
		  overflow: auto; /* Enable scroll if needed */
		  background-color: rgb(0,0,0); /* Fallback color */
		  background-color: rgba(0,0,0,0.4); /* Black w/ opacity */
		}

		/* Modal Content/Box */
		.modal-content {
		   position:fixed;
		  background-color: #fefefe;
		  margin: 15% auto; /* 15% from the top and centered */
		  padding: 20px;
		  left:10%;
		  border: 1px solid #888;
		  width: 80%; /* Could be more or less, depending on screen size */
		  word-break: break-all;
		  z-index:999;
		}

		/* The Close Button */
		.close {
		  color: #aaa;
		  float: right;
		  font-size: 28px;
		  font-weight: bold;
		}
		hr { 
		    display: block;
		    margin-before: 0.5em;
		    margin-after: 0.5em;
		    margin-start: auto;
		    margin-end: auto;
		    overflow: hidden;
		    border-style: inset;
		    border-width: 1px;
		}
		a.tooltip {
	        user-select: none;
	        -moz-user-select: none;
	        -khtml-user-select: none;
	        -webkit-user-select: none;
	        -o-user-select: none;
		  position: relative;
		  text-decoration: none;
		}
.tooltip {
	position:fixed;
	        user-select: none;
	        -moz-user-select: none;
	        -khtml-user-select: none;
	        -webkit-user-select: none;
	        -o-user-select: none;
			  position: absolute;
			  padding:10px;
			    right:20px;
			  top:0px;
			  height:5px;
			  font-size:10px;
			  z-level:0;
			  background-color:yellow;
			  border-radius:0px 0px 15px 15px;
			  animation-delay: 10s;
			    -webkit-animation: fadein 2s; /* Safari, Chrome and Opera > 12.1 */
			       -moz-animation: fadein 2s; /* Firefox < 16 */
			        -ms-animation: fadein 2s; /* Internet Explorer */
			         -o-animation: fadein 2s; /* Opera < 12.1 */
			            animation: fadein 2s;
}

@keyframes fadein {
    from { opacity: 0; }
    to   { opacity: 1; }
}

/* Firefox < 16 */
@-moz-keyframes fadein {
    from { opacity: 0; }
    to   { opacity: 1; }
}

/* Safari, Chrome and Opera > 12.1 */
@-webkit-keyframes fadein {
    from { opacity: 0; }
    to   { opacity: 1; }
}

/* Internet Explorer */
@-ms-keyframes fadein {
    from { opacity: 0; }
    to   { opacity: 1; }
}

/* Opera < 12.1 */
@-o-keyframes fadein {
    from { opacity: 0; }
    to   { opacity: 1; }
}

</style>
  </head>
  <body>';

#
## Title
#

print '<h3>Proteinortho 6 : '.$proteinorthotsvfile.'</h3>';

print '<div id="myModal" class="modal">

  <!-- Modal content -->
  <div class="modal-content">
    <span class="close">&times;</span>
    <h3 style="color:darkblue;">Extract fasta from selected group</h3>
    <hr></hr>
    <p id="modaltxt">Some text in the Modal..</p>
  </div>

</div>
<div class="tooltip" >Tip: You can click on the first 3 columns of a row (blue) to get information how to extract the genes/proteins from the given fastas.</div>
';
#
## Navigation Bar 
#

print '<div class="overlay" id="overlay"><div class="left"><input type="text" id="searchField" onchange="isEmpty()" onkeyup="this.onchange();" onpaste="this.onchange();" oninput="this.onchange();" placeholder="Search..">';
print '<br><input type="button" onclick="filterTable(0)" value="search" id="searchbutton">';#<input type="button" onclick="filterTable(1)" value="exact search">
print '<br><span>The search for <span id="countspanA"></span><br>results in <b><span id="countspanB"></span></b> groups.</span></div><div class="middleA"><div class="sidewaysdiv">manage</div><input type="button" onclick="toggleCols();" value="toggle columns" id="toggleColsButton"><br><input type="button" onclick="togglePlus();" value="toggle (+)" id="togglePlusButton"><br><input type="checkbox" id="toggleShowTrueNames" name="toggleShowTrueNames"
          onclick="toggleShowTrueNames();"><label id="toggleShowTrueNamesLabel" for="toggleShowTrueNames">show IDs only</label></div>
<div class="middleB"><div class="sidewaysdiv">select</div><select multiple data-placeholder="Begin typing a name to filter..." class="chosen-select no-scroll" name="test" id="select">
';
	open(FILE,"<",$proteinorthotsvfile);
	while(<FILE>){
		my $curline=$_;
		my $i=0;
		foreach $k (split("\t",$curline)){
			if($i<3){
			}else{
				if($k=~m/([^_]+)_([^_]+)_(.*)/){
					$l=$1;$m=$2;$n=$3;
					$n=~s/.f[an]a//;
					$k=$l."_".$m." ".$n;
				}
				print "<option selected>&#x2611; $k</option>";
			}
			$i++;
		}
		last;
	}
	close(FILE);
print '</select></div></div>';

#
## TABLE
#

print '<div class="clusterize"><div id="scrollArea" class="clusterize-scroll"><table id="myTable" class="greyGridTable"><thead><tr class="data-heading">';

	open(FILE,"<",$proteinorthotsvfile);
	while(<FILE>){
		my $curline=$_;
		my $i=0;
		foreach $k (split("\t",$curline)){
			if($i<3){
				print "<th>$k</th>";# onclick='sortTable($i)'
			}else{
				if($k=~m/([^_]+)_([^_]+)_(.*)/){
					$l=$1;$m=$2;$n=$3;
					$n=~s/.f[an]a//;
					$k=$l."_".$m." ".$n;
				}
				print "<th>$k</th>";
			}
			$i++;
		}
		last;
	}
	close(FILE);
	print '</tr></thead>';


	print '<tbody id="contentArea" class="clusterize-content">
        <tr class="clusterize-no-data">
          <td>Loading data…</td>
        </tr>
      </tbody>';


	print '<tfoot><tr>';

	open(FILE,"<",$proteinorthotsvfile);
	while(<FILE>){
		my $curline=$_;
		my $i=0;
		foreach $k (split("\t",$curline)){
			if($i<3){
				print "<th>$k</th>";
			}else{
				if($k=~m/([^_]+)_([^_]+)/){
					$k=substr($1,0,1).".".substr($2,0,3);
				}else{
					$k=substr($k,0,5);
				}
				print "<th>".$k."</th>";
			}
			$i++;
		}
		last;
	}
	close(FILE);
	print '</tr>';
	print '</tfoot>';

print '</table></div></div>';


print '</body>';

#
## JS part2
#

print '<script type="text/javascript">
		var input = document.getElementById("searchField");
		// Execute a function when the user releases a key on the keyboard
		input.addEventListener("keyup", function(event) {
		  // Number 13 is the "Enter" key on the keyboard
		  if (event.keyCode === 13) {
		    // Cancel the default action, if needed
		    event.preventDefault();
		    // Trigger the button element with a click
		    document.getElementById("searchbutton").click();
		  }
		}); 

		var selectedCols=[];
		function getselectValues(select) {
		  var result = [];
		  var options = select && select.options;
		  var opt;
	      result.push(0);
		  result.push(1);
		  result.push(2);
		  
		  for (var i=0, iLen=options.length; i<iLen; i++) {
		    opt = options[i];

		    if (opt.outerHTML.includes("selected")) {
		      result.push(i+3);
		    }
		  }
		  return result;
		}

		select=document.getElementById(\'select\');
		select.addEventListener("mouseenter", mouseOver);
		select.addEventListener("mouseleave", mouseOut);
		var B = document.body,
		    H = document.documentElement,
		    height

		if (typeof document.height !== \'undefined\') {
		    height = document.height // For webkit browsers
		} else {
		    height = Math.max( B.scrollHeight, B.offsetHeight,H.clientHeight, H.scrollHeight, H.offsetHeight );
		}
		function mouseOver() {
			// table=document.getElementById(\'myTable\');
			//table.style.zIndex="-1";
			//select.style.zIndex="999";
			//select.style.height=(height-50)+"px";
			//overlay=document.getElementById(\'overlay\');
			//overlay.style.zIndex="1";

		}

		function mouseOut() {
			//table=document.getElementById(\'myTable\');
			//table.style.zIndex="999";
			//select.style.zIndex="1";
			//select.style.height="";
			//overlay=document.getElementById(\'overlay\');
			//overlay.style.zIndex="-1";

		}


		window.onmousedown = function (e) {
		    var el = e.target;
		    if (el.tagName.toLowerCase() == \'option\' && el.parentNode.hasAttribute(\'multiple\')) {
		        e.preventDefault();

		        // toggle selection
		        if (el.hasAttribute(\'selected\')){ 
		        	el.removeAttribute(\'selected\');
		        	var text = el.textContent || el.innerText;
		        	el.innerHTML = \'&#x2610; \'+text.substring(1);
		        }else{ 
		        	el.setAttribute(\'selected\', \'\');
		        	var text = el.textContent || el.innerText;
		        	el.innerHTML = \'&#x2611; \'+text.substring(1);
		        }

		        // hack to correct buggy behavior
		        //var select = el.parentNode.cloneNode(true);
		        //el.parentNode.parentNode.replaceChild(select, el.parentNode);

		        updateTable();
		    }
		}

		function toggleCols() {
			var selects=document.getElementsByTagName(\'select\')[0];
		 	var optionsV = selects.options;
			
			for (i = 0; i < optionsV.length; i++) {
		        if (optionsV[i].hasAttribute(\'selected\')){ 
		        	optionsV[i].removeAttribute(\'selected\');
		        	var text = optionsV[i].textContent || optionsV[i].innerText;
		        	optionsV[i].innerHTML = \'&#x2610; \'+text.substring(1);
		        }else{ 
		        	optionsV[i].setAttribute(\'selected\', \'\');
		        	var text = optionsV[i].textContent || optionsV[i].innerText;
		        	optionsV[i].innerHTML = \'&#x2611; \'+text.substring(1);
		        }
		    }

			updateTable();
		}

		function togglePlus() {
			var selects=document.getElementsByClassName(\'popup\');
			for (i = 0; i < selects.length; i++) {
				selects[i].click();
		    }
		}

		var data_raw = [';

		open(FILE,"<",$proteinorthotsvfile);
		my $rowi=0;
		while(<FILE>){
			if($rowi==0){$rowi++;next;}
			if($rowi>1){print ',';}
			chomp;
			my $curline=$_;
			print '\'<tr>';
			my $coli=0;
			foreach $k (split("\t",$curline)){
				if($coli>2){ 
					my @ar=split(",",$k); 
					
					print "<td>".$ar[0];
					
					if(scalar(@ar)>1){

						$hiddentext="";
						for(my $ii=1;$ii<scalar(@ar);$ii++){
							if($hiddentext ne ""){
								$hiddentext.=", ";
							}
							$hiddentext.=$ar[$ii];
						}
						print '<div class="popup" onclick="showpopup(\\\''.$rowi."_".$coli.'\\\')" id="myPopup'.$rowi."_".$coli.'"><font color="blue">(+)</font></div><div class="popuptext" id="myPopupText'.$rowi."_".$coli.'"> '.$hiddentext.'</div>';
					}
					print "</td>"; 
				}else{
					if($coli < 3){
						print "<td onclick=\"tableText(this.parentNode)\">$k</td>";
					}else{
						print "<td>$k</td>";
					}
				}
				$coli++;
			}
			print '</tr>\'';
			$rowi++;
		}
		close(FILE);
print '];

		var data_raw_mapped = [';

		if(scalar(keys(%mappingdata))>0){
			open(FILE,"<",$proteinorthotsvfile);
			my $rowi=0;
			while(<FILE>){
				if($rowi==0){$rowi++;next;}
				if($rowi>1){print ',';}
				chomp;
				my $curline=$_;
				print '\'<tr>';
				my $coli=0;
				foreach $k (split("\t",$curline)){
					if($coli>2){ 
						my @ar=split(",",$k);

						if(exists $mappingdata{$ar[0]}){
							$tmp=$mappingdata{$ar[0]};
							$tmp=~s/[^a-zA-Z0-9. \-]//gm;
							print "<td>".$tmp."(".$ar[0].")";
						}else{
							print "<td>".$ar[0];
						}

						if(scalar(@ar)>1){

							$hiddentext="";
							for(my $ii=1;$ii<scalar(@ar);$ii++){
								if($hiddentext ne ""){
									$hiddentext.=", ";
								}
								if(exists $mappingdata{$ar[$ii]}){
									$tmp=$mappingdata{$ar[$ii]};
									$tmp=~s/[^a-zA-Z0-9. \-]//gm;
									$hiddentext.=$tmp."(".$ar[$ii].")";
								}else{
									$hiddentext.=$ar[$ii];
								}
							}
							print '<div class="popup" onclick="showpopup(\\\''.$rowi."_".$coli.'\\\')" id="myPopup'.$rowi."_".$coli.'"><font color="blue">(+)</font></div><div class="popuptext" id="myPopupText'.$rowi."_".$coli.'"> '.$hiddentext.'</div>';
						}
						print "</td>"; 
					}else{
						if($coli < 3){
							print "<td onclick=\"tableText(this.parentNode)\">$k</td>";
						}else{
							print "<td>$k</td>";
						}
					}
					$coli++;
				}
				print '</tr>\'';
				$rowi++;
			}
			close(FILE);
		}
		
print '];

		var useMappedData=false;

		if(data_raw_mapped.length==0){
			document.getElementById(\'toggleShowTrueNames\').style.display="none";
			document.getElementById(\'toggleShowTrueNamesLabel\').style.display="none";
		}else{
			useMappedData=true;
		}

		function toggleShowTrueNames() {
			useMappedData=!useMappedData;
			updateTable(true);
		}

		// sort by number of species:
			var data_first_then_second_then_third_col_sort_idx=[];
			for (j = 0; j < data_raw.length; j++) {
				data_first_then_second_then_third_col_sort_idx[j]=j;
			}
			data_first_then_second_then_third_col_sort_idx=data_first_then_second_then_third_col_sort_idx.sort(function(a, b) {

			  aa=(data_raw[a].replace(/^<tr><td[^>]*>/,"").split(/<\/td><td[^>]*>/))[0];
			  bb=(data_raw[b].replace(/^<tr><td[^>]*>/,"").split(/<\/td><td[^>]*>/))[0];
			  if(aa != bb){
			  	return bb - aa;
			  }else{
			  	aa=(data_raw[a].replace(/^<tr><td[^>]*>/,"").split(/<\/td><td[^>]*>/))[1];
			  	bb=(data_raw[b].replace(/^<tr><td[^>]*>/,"").split(/<\/td><td[^>]*>/))[1];
			  	if(aa != bb){
				  	return bb - aa;
				  }else{
				  	aa=(data_raw[a].replace(/^<tr><td[^>]*>/,"").split(/<\/td><td[^>]*>/))[2];
				  	bb=(data_raw[b].replace(/^<tr><td[^>]*>/,"").split(/<\/td><td[^>]*>/))[2];
				  	return bb - aa;
				  }
			  }

			});
			data_raw = data_first_then_second_then_third_col_sort_idx.map(i => data_raw[i]);
			if(data_raw_mapped.length>0){
				var data_first_then_second_then_third_col_sort_idx=[];
				for (j = 0; j < data_raw_mapped.length; j++) {
					data_first_then_second_then_third_col_sort_idx[j]=j;
				}
				data_first_then_second_then_third_col_sort_idx=data_first_then_second_then_third_col_sort_idx.sort(function(a, b) {

				  aa=(data_raw_mapped[a].replace(/^<tr><td[^>]*>/,"").split(/<\/td><td[^>]*>/))[0];
				  bb=(data_raw_mapped[b].replace(/^<tr><td[^>]*>/,"").split(/<\/td><td[^>]*>/))[0];
				  if(aa != bb){
				  	return bb - aa;
				  }else{
				  	aa=(data_raw_mapped[a].replace(/^<tr><td[^>]*>/,"").split(/<\/td><td[^>]*>/))[1];
				  	bb=(data_raw_mapped[b].replace(/^<tr><td[^>]*>/,"").split(/<\/td><td[^>]*>/))[1];
				  	if(aa != bb){
					  	return bb - aa;
					  }else{
					  	aa=(data_raw_mapped[a].replace(/^<tr><td[^>]*>/,"").split(/<\/td><td[^>]*>/))[2];
					  	bb=(data_raw_mapped[b].replace(/^<tr><td[^>]*>/,"").split(/<\/td><td[^>]*>/))[2];
					  	return bb - aa;
					  }
				  }

				});
				data_raw_mapped = data_first_then_second_then_third_col_sort_idx.map(i => data_raw_mapped[i]);
			}

		function updateTable(re_draw=false){

			selectedCols=getselectValues(document.getElementsByTagName(\'select\')[0]);
			
			table = document.getElementById("myTable");
			tr = table.getElementsByTagName("tr");
			tds = table.getElementsByTagName("td");
			num_cols=tds.length/(tr.length-3);
			num_rows=tr.length;

			if(re_draw || tr.length<4){ //empty clusterize table -> init

				if(useMappedData){
					if(!clusterize){
						clusterize = new Clusterize({
						  rows: data_raw_mapped,
						  scrollId: \'scrollArea\',
						  contentId: \'contentArea\'
						});
					}else{
						clusterize.update(data_raw_mapped);
					}
				}else{
					if(!clusterize){
						clusterize = new Clusterize({
						  rows: data_raw,
						  scrollId: \'scrollArea\',
						  contentId: \'contentArea\'
						});
					}else{
						clusterize.update(data_raw);
					}
				}
				
				table = document.getElementById("myTable");
				tr = table.getElementsByTagName("tr");
				tds = table.getElementsByTagName("td");
				num_cols=tds.length/(tr.length-3);
				num_rows=tr.length;

			}else{
				var resultArr;
				if(useMappedData){
					resultArr = data_filtered_idx.map(i => {
						dat=data_raw_mapped[i].substring(0,data_raw_mapped[i].length-10).split(/<\/td><td[^>]*>/); //length-10 => "</td></tr>" in the last column.
						return selectedCols.slice(0,3).map(j => dat[j]).join("</td><td onclick=\"tableText(this.parentNode)\">") + "</td><td>" + selectedCols.slice(3).map(j => dat[j]).join("</td><td>")+"</td></tr>";
					}); 
				}else{
					resultArr = data_filtered_idx.map(i => {
						dat=data_raw[i].substring(0,data_raw[i].length-10).split(/<\/td><td[^>]*>/); //length-10 => "</td></tr>" in the last column.
						return selectedCols.slice(0,3).map(j => dat[j]).join("</td><td onclick=\"tableText(this.parentNode)\">") + "</td><td>" + selectedCols.slice(3).map(j => dat[j]).join("</td><td>")+"</td></tr>";
					}); 
				}
				clusterize.update(resultArr);
			}

			if(tr.length>3){
				var tds=tr[0].getElementsByTagName("th");
				for (j = 0; j < tds.length; j++) {
					if( tds[j] ){
						if(selectedCols.indexOf(j)==-1){
							(tds[j]).style.display="none";
						}else{
							(tds[j]).style.display="";
						}
					}
				}
				var tds=tr[tr.length-1].getElementsByTagName("th");
				for (j = 0; j < tds.length; j++) {
					if( tds[j] ){
						if(selectedCols.indexOf(j)==-1){
							(tds[j]).style.display="none";
						}else{
							(tds[j]).style.display="";
						}
					}
				}
			}

		}

		isEmpty();

		var table = document.getElementById("myTable");

';
if(scalar(@ARGV)>1){
	print 'var fastas ="';
	# test if there is a common prefix
	my $foundCommonPrefix=0;
	my $commonPrefix="";
	my $commonSuffix="";
	if(scalar(@ARGV)>1){
		$foundCommonPrefix=1;
		$commonPrefix=$proteinorthotsvfile;
		$commonSuffix=$proteinorthotsvfile;
		$commonPrefix=~s/\/[^\/]*$/\//;
		$commonSuffix=~s/.*\.([^.]+)$/$1/;
		for(my $argi=2;$argi<scalar(@ARGV);$argi++){
			my $curPrefix=$ARGV[$argi];
			my $curSuffix=$ARGV[$argi];
			$curPrefix=~s/\/[^\/]*$/\//;
			$curSuffix=~s/.*\.([^.]+)$/$1/;
			if($commonPrefix ne $curPrefix || $commonSuffix ne $curSuffix){
				$foundCommonPrefix=0;
				break;
			}
		}	
	}
	#if not print all
	if(!$foundCommonPrefix){
		for(my $argi=1;$argi<scalar(@ARGV);$argi++){
			if($argi!=1){
				print ' ';
			}
			print $ARGV[$argi];
		}	
	}else{
		print $commonPrefix."*".$commonSuffix;
	}
	print '";';
}else{
	print 'var fastas ="INPUT_FASTAS";';
}
print '
		var modal = document.getElementById(\'myModal\');
		var modaltxt = document.getElementById(\'modaltxt\');
		var span = document.getElementsByClassName("close")[0];
		window.onclick = function(event) {
		  if (event.target == modal) {
		    modal.style.display = "none";
		  }
		} 
		span.onclick = function() {
		  modal.style.display = "none";
		}
		function tableText(tableCell) {
			var tds=tableCell.getElementsByTagName("td");
			var genes="";
		    for (var i = 3; i < tds.length; i++) {
		    	if(tds[i].innerText!="*"&&tds[i].innerText!=""){
		    		if(genes!="" && genes.slice(-1)!=","){
		    			genes+=",";
		    		}
		    		genes+=tds[i].innerHTML.replace(/([<]([^>]+)[>])/ig,"").replace(/[^()]+[(]([^)]+)[)]/ig,"$1,").replace(/ /ig,"").replace(/,$/ig,"").replace(/([(]([+-]+)[)])/ig,",").replace(/,,+/ig,",").replace(/[()]/g,"");
		    	}
		    }

		    modal.style.display = "block";
		    modaltxt.innerHTML="Use the following command to obtain the fasta-file of the selected group:<br><span><br>proteinortho_grab_proteins.pl -exact \'"+genes+"\' "+fastas+"</span><br><br>";
		    if(fastas=="INPUT_FASTAS"){
		    	modaltxt.innerHTML+="Replace INPUT_FASTAS with input fasta files (e.g. test/*faa).<br>";
		    }
		    modaltxt.innerHTML+="Use the appendix \'>SOMEFILE.fasta\' to redirect the output to SOMEFILE.fasta.<br>You can use the option \'-source\' for adding the species name to each gene in the output files.<br>(If you are missing \'poteinortho_grab_proteins.pl\' you can download it <a href=\"https://gitlab.com/paulklemm_PHD/proteinortho/blob/master/src/proteinortho_grab_proteins.pl\">here</a>)";
		}
		
		</script></html>';

