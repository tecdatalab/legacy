<?php
include ('../header.php');
?>
<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/2.6.0/Chart.min.js">
    </script>
<br>

<?php

$db_prog = 'query_top_n.pl';
$img_dir = '../db/img';
$db_dir = '../db';

$search_upload = 0;
$query_id = $_REQUEST['q'];
$volume_filter = $_REQUEST['volumeFilter'];
$min_resolution = $_REQUEST['minResolutionFilter'];
$max_resolution = $_REQUEST['maxResolutionFilter'];
$representation = $_REQUEST['rep'];

#$ccc = count($_REQUEST);
#echo "Number $ccc<br>";

#foreach (array_keys($_REQUEST) as $currkey) {
#	echo "$currkey <br>";
#}
// Negative values are assumed to be ignored
if (!$min_resolution || !($min_resolution >= 0)) {
	$min_resolution = -1;
}
if (!$max_resolution || !($max_resolution >= 0)) {
	$max_resolution = -1;
}

$inv_dir = '../db/2_contour/inv';
if ($representation === "recommend") {
	$inv_dir = '../db/recommend_contour/inv';
}
if ($representation === "recommendonethird") {
	$inv_dir = '../db/1_contour/inv';
}
if ($representation === "recommendtwothird") {
	$inv_dir = '../db/2_contour/inv';
}
if ($representation === "recommendonethirdtwothird") {
	$inv_dir = '../db/1_2_contour/inv';
}
if ($representation === "recommendonestd") {
	$inv_dir = '../db/1_std_merge/inv';
}

if ($volume_filter === "off") {
	$volume_filtering_option = "-nofilter";
} else {
	$volume_filtering_option = "-filter";
}

# Load all hashtables from configuration files
$P2C_file = fopen("$db_dir/name.txt", "r");
while (!feof($P2C_file)) {
	$line = fgets($P2C_file);
	#separate by a char
	$tmparray = preg_split("/\t/", $line);
	$namehash["$tmparray[0]"] = "$tmparray[1]";
}
fclose($P2C_file);

$P2C_file = fopen("$db_dir/fullname.txt", "r");
while (!feof($P2C_file)) {
	$line = fgets($P2C_file);
	#separate by a char
	$tmparray = preg_split("/\t/", $line);
	$fullnamehash["$tmparray[0]"] = "$tmparray[1]";
}
fclose($P2C_file);

$P2C_file = fopen("$db_dir/volume.txt", "r");
while (!feof($P2C_file)) {
	$line = fgets($P2C_file);
	#separate by a char
	$tmparray = preg_split("/\t/", $line);
	$volumehash["$tmparray[0]"] = "$tmparray[1]";
}
fclose($P2C_file);

$P2C_file = fopen("$db_dir/resolutions.txt", "r");
while (!feof($P2C_file)) {
	$line = fgets($P2C_file);
	#separate by a char
	$tmparray = preg_split("/\t/", $line);
	$resolutionhash["$tmparray[0]"] = "$tmparray[1]";
}
fclose($P2C_file);

// Configure INV files and names depending on upload or regular
// EMDB ID search
$uploaded_tmp_file = $_FILES['pdbfile']['tmp_name'];
if ($uploaded_tmp_file) {
	$search_upload = 1;
	// Original name of the file submitted
	$query_name = $query_fullname = $_FILES['pdbfile']['name'];
	$query_file = "$uploaded_tmp_file.inv";
	$contour_level = $_REQUEST['contour'];
#	$upload_log = `uploadscripts/em_to_3dzd.sh $uploaded_tmp_file $contour_level $query_file`;
	$upload_log = `uploadscripts/upload_to_inv.pl $uploaded_tmp_file $contour_level $representation $query_file`;
	// Compute the query_volume. The second position in the text returned is the volume.
	$volume_result = explode(" ", `uploadscripts/em_volume $contour_level $uploaded_tmp_file`);
	$query_volume = $volume_result[1];
} else { // Search using an existing EMDB ID
	$query_name = $namehash["$query_id"];
	$query_fullname = $fullnamehash["$query_id"];
	$query_volume = $volumehash["$query_id"];
	$query_file = "$inv_dir/$query_id.inv";
}

// Generate the text results URL
$text_results_url = "listResults.cgi";
if ($uploaded_tmp_file) {
	$formatted_query_file = urlencode(str_replace("/tmp/", "", $query_file));
	$text_results_url .= "?mapfile=$formatted_query_file&mapvolume=" . urlencode($query_volume);
} else {
	$text_results_url .= "?emdbid=$query_id";
}
$text_results_url .= "&volumefilter=$volume_filter&representation=$representation";
if ($min_resolution >= 0) {
	$text_results_url .= "&minresolution=" . urlencode($min_resolution);
}
if ($max_resolution >= 0) {
	$text_results_url .= "&maxresolution=" . urlencode($max_resolution);
}
// show the inv curve and values
$zernike_length_end_pos;
$zernike_descriptors;
$zernike_count;

// ---- Check INV File ----
$cat_inv = file_get_contents($query_file);

// initialize the zernike descriptors variable that will be used at the end of the table
$zernike_length_end_pos = strpos($cat_inv, "\n");
$zernike_descriptors = substr($cat_inv, $zernike_length_end_pos);

$zernike_count = 0;
$descriptor;
$chart_formatted_descriptors = "";
$textarea_formatted_descriptors = "";

$tmparray = preg_split("/\n/", $zernike_descriptors);
for ($i = 0; $i < count($tmparray) - 1; ++ $i) {
	$descriptor = $tmparray[$i];
	if ($zernike_count > 0) {
		$textarea_formatted_descriptors = "$textarea_formatted_descriptors($zernike_count) $descriptor \n";
		$chart_formatted_descriptors = "$chart_formatted_descriptors$zernike_count,$descriptor;";
	}
	$zernike_count++;
}

// All error messages should be handled here. There will be a PHP IF clause
// whose scope is the rest of the HTML code generated for visualization.
// If any error case is found then just the error message should be displayed,
// followed by the footer page.
// Thus there will be one closing } at the very end of this php file
// to close out this if-else error check
include ('errortable.php');
// Error #1: The EMDB ID provided does not exist overall
if (!$search_upload && !file_exists($query_file) && $query_fullname == "")
{
	generate_error_table("Error: Entry not in current EMDB.", "The EMDB ID provided is not included in the current version of the database.");
include ('../tailer.php');
}
// Error #2: The EMDB ID exists but there is some error in the generation of 3DZDs
else if(!$search_upload && !file_exists($query_file))
{
	generate_error_table("Error: EMDB ID in process", "The EMDB ID provided is part of the current update but all necessary files are not generated yet.<br>" .
							  "A fraction of EMDB entries present problems either because they require manual formatting to be used by our " .
							  " programs or, in some cases, the data in EMDB servers is incomplete.");
include ('../tailer.php');
}
else
{
?>


<!-- summarize the query -->
<table width="768" align="center" border="1">
	<tr>
		<td colspan="4"><div align="left"><p class="text4">Query entry
		(<a href=<?php echo "$text_results_url"?>>Download text results here</a>)</p></div></td>
	</tr>
	<tr>
		<td width="192">
		<div align="left">
<!-- When it's an uploaded file, the left panel only shows the name of the file -->
<?php if ($search_upload) { ?>
		<p class="text4">User uploaded query</p>
		<p><b>File uploaded:</b></p>
		<p><?php echo $query_name ?></p>
<!-- When it's a regular search, it displays a link, image and description -->
<?php } else { ?>
		<b>
		<a
			href='http://www.ebi.ac.uk/pdbe/entry/EMD-<?php echo "$query_id" ?>'
			target="_blank" title='<?php echo "EMD-$query_id" ?>'>
			<?php echo "EMD-$query_id" ?>
		</a>
		</b>
		<br>
		<?php if ($query_fullname != $query_name) { ?>
		<a title='<?php echo "$query_fullname" ?>'> 
		<?php } ?>
		<font color="black"><?php echo "$query_name" ?></font>
		<?php if ($query_fullname != $query_name) { ?>
		</a>
		<?php } ?>
		<br>
		
		<img alt='<?php echo "$query_id" ?>' width="180" border="0" src='<?php echo "$img_dir/$query_id.gif" ?>' />	
		<br>
<?php } ?>
		</div>
		</td>
		
		<!-- create the applet that displays the zernike descriptors -->
		<td colspan="2" width="384">
    		
<!--
    		<applet name="chart" code="JFreeChartApplet.class"
				codebase="../jfreechart"
				archive="jcommon-1.0.15.jar,jfreechart-1.0.12.jar" width="384">
				<PARAM name="points"
					value='<?php echo "$chart_formatted_descriptors"; ?>' />
			</applet>
-->
            <div style="max-width: 400px; max-height: 200px">
                <canvas id="chart-js" width="400" height="200"></canvas>
            </div>
		</td>
			
		<td width="192">
		    <b>Descriptors</b>
			<br>
			<textarea align = "left" id="statusbarzernike" rows="11" readonly><?php echo "$textarea_formatted_descriptors"; ?>
			</textarea>
		</td>	
			
	</tr>
	
	<tr>
		<td colspan="4"><div align="left"><p class="text4">Retrieval results with similar 3D surface
		(<?php echo "Volume Filter: $volume_filter";
		if ($min_resolution >= 0) echo ", Min. Resolution: $min_resolution";
		if ($max_resolution >= 0) echo ", Max. Resolution: $max_resolution";
		?>)
		</p></div></td>
	</tr>

		<!-- Render the Surface Retrieval Results -->
<?php
if($search_upload) {
	$query_result = `./$db_prog 20 $query_file $volume_filtering_option $inv_dir $min_resolution $max_resolution -upload $query_volume`;
} else {
	$query_result = `./$db_prog 20 $query_id $volume_filtering_option $inv_dir $min_resolution $max_resolution`;
}
$tmparray = preg_split("/\n/", $query_result);
?>

<tr>
		<?php

$count = 0;

for ($i = 0; $i < count($tmparray) - 1; ++ $i) {

	++ $count;
	$resultarray = preg_split("/\t/", $tmparray[$i]);
	$template_id = $resultarray[0];
	$dist = $resultarray[1];
	$template_name = $namehash["$template_id"];
	$template_fullname = $fullnamehash["$template_id"];
	$template_volume = $volumehash["$template_id"];
	$template_resolution = $resolutionhash["$template_id"];
	if (!$template_resolution) {
		$template_resolution = "N/A";
	}
	if (($query_volume > 0) && ($template_volume > 0)) {
		$ratio_volume = $template_volume / $query_volume;
		$ratio_volume = sprintf("%.3f", $ratio_volume);
	} else {
		$ratio_volume = "N/A";
	}
	$imghtmlink = "./search.php?q=" . $template_id . "&volumeFilter=" . $volume_filter . "&rep=" . $representation;
?>
	
		
		<td width="192">
		<b>
		<a
			href='http://www.ebi.ac.uk/pdbe/entry/EMD-<?php echo "$template_id" ?>'
			target="_blank" title='<?php echo "EMD-$template_id" ?>'>
			<?php echo "EMD-$template_id" ?>
		</a>
		</b>
		<br>
		<?php if ($template_fullname != $template_name) { ?>
		<a title='<?php echo "$template_fullname" ?>'> 
		<?php } ?>
		<font color="black"><?php echo "$template_name" ?></font>
		<?php if ($template_fullname != $template_name) { ?>
		</a>
		<?php } ?>
		<br>
		<a class="imgLink" href='<?php echo "$imghtmlink" ?>'>
		<img alt='<?php echo "$template_id" ?>' width="180" border="0" src='<?php echo "$img_dir/$template_id.gif" ?>' />
		</a>
		<br>
		<font color="blue"><?php echo "EucD: $dist" ?></font>
		<br>
		<font color="blue"><?php echo "Ratio of volume: $ratio_volume" ?></font>
		<br>		
		<font color="blue"><?php echo "Resolution: $template_resolution" ?></font>
		<br>		
		</td>
		
		<?php


	if (($count % 4 == 0) && $count != 20) {
?>
	</tr>
	<tr>
	<?php


	} else
		if ($count == 20) {
?>
</tr>
<?php
		}
}

?>

</table>
<script>
//for chart
var ctx = document.getElementById("chart-js").getContext('2d');
console.log(ctx);
var descriptors = "<?php echo "$chart_formatted_descriptors"; ?>";
descriptors = descriptors.split(";");
var labels = [];
for (var i = 0; i < descriptors.length; i ++) {
    labels[i] = descriptors[i].split(",")[0];
    descriptors[i] = descriptors[i].split(",")[1];
}

var chart = new Chart(ctx, {
    // The type of chart we want to create
    type: 'line',

    // The data for our dataset
    data: {
        labels: labels,
        datasets: [{
            label:"Zernike",
            backgroundColor: 'rgba(255, 255, 255, 0)',
            borderColor: 'rgb(20, 20, 20)',
            borderWidth: 1,
            data: descriptors
        }],
    },
    scaleLineColor: 'red',
    animation: true,
    // Configuration options go here
    options : {
        animation:{
            onComplete : function(){
                try{

                } catch(e) {
                    
                }
            }
        },
        legend: {
            display: false
        },
        
        elements: { point: { radius: 0 } },
        responsive: true,
        maintainAspectRatio: false,
        scales: {
            yAxes: [{
              scaleLabel: {
                display: true,
                labelString: 'Value',
                fontSize: 20,
                fontColor: '#000'
              },
              ticks: {
                fontColor: '#000',
                fontSize: 15
              }
              
            }],
            xAxes: [{
              scaleLabel: {
                display: true,
                labelString: 'Zernike Descriptor Number',
                fontSize: 20,
                fontColor: '#000'
              },
              ticks: {
                autoSkip: true,
                maxTicksLimit: 10,
                maxRotation: 0,
                fontColor: '#000',
                fontSize: 15
              }
            }]
          }
        }
});

  
</script>
<?php


include ('../tailer.php');
} // CLOSE normal visualization generation (coming from the error management part)
?>
