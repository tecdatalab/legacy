<?php
include ('../header.php');
?>

<?php

$query_pdb_id_list = $_POST['listarea'];
$query_pdb_file = $_FILES['listfile']['tmp_name'];
$volume_filter = $_POST['volumeFilter'];
$num = $_POST['num'];
$representation = $_REQUEST['rep'];

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

if ($volume_filter === "off") {
	$volume_filtering_option = "-nofilter";
} else {
	$volume_filtering_option = "-filter";
}

if ($query_pdb_file) {

	$unique_filename = substr($query_pdb_file, 5);
	#move
	$target = "tmpbatch/" . $unique_filename;
	move_uploaded_file($query_pdb_file, $target);

	# read the content of the uploaded file
	$fh = fopen($target, 'r');
	$thefiledata = fread($fh, filesize($target));
	fclose($fh);

	$query_pdb_id_list = $thefiledata;
}

$pattern = "/:|\s+|\n/";
$query_pdb_id_list = preg_replace($pattern, ",", $query_pdb_id_list);

$randnum = rand();
$time = date("Him");
$queryresult_file_path = "queryresult/".$randnum . $time;
$queryresult_file_path_zip = $queryresult_file_path . ".zip";

`./searchbatch.pl $num $volume_filtering_option $query_pdb_id_list $inv_dir $queryresult_file_path`;
`zip -r $queryresult_file_path_zip $queryresult_file_path`;
?>

<script
	type='text/javascript' src='../jsFunctions.js'></script>
<script
	type='text/javascript' src='../ajaxFunctions.js'></script>
<script
	type='text/javascript' src='../jquery.all.js'></script>
<link
	rel="stylesheet" type="text/css" href="../jquery.autocomplete.css" />

<table width="100%" border="0" cellspacing="0" cellpadding="0">
	<tr>
		<td>
			<table width="700" border="0" align="center" cellpadding="0"
				cellspacing="0">
				<tr valign="top">
					<td width="700">
						<table width="100%" border="0" cellspacing="0" cellpadding="0">
							<tr>
								<td><img src="../images/spacer.gif" alt="" width="100"
									height="35" />
								</td>
							</tr>
							<tr>
								<td class="text4">Result shown</td>
							</tr>
							<tr>
								<td><p>&nbsp;</p>
								<p><a href="<?php echo $queryresult_file_path_zip; ?>"><H2>Please download all result files in a compressed package</H2></a></p>
								<br>
								<!-- <p><?php echo $query_pdb_id_list; ?></p> -->
								
								 <?php

$d = dir($queryresult_file_path) or die("Wrong path: $queryresult_file_path");
while (false !== ($entry = $d->read())) {
	if ($entry != '.hit' && $entry != '.' && $entry != '..' && !is_dir($dir . $entry)) {
		$filepath = $queryresult_file_path. "/".$entry;
?>
                                 <p><a href="<?php echo $filepath; ?>"><?php echo $entry; ?></a></p>
                           
                                 <?php
	}
}
$d->close();
?>                                
                                 <p>&nbsp;</p>
                                 <p>&nbsp;</p>
                                 <p>&nbsp;</p>
                                </td>
							</tr>
						</table>
					</td>
				</tr>
			</table>
		</td>
	</tr>
</table>


<!-- surmarize the query -->
<?php


include ('../tailer.php');
?>
