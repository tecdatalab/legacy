<?php
include ('header.php');
?>

<?php

$latest_updated_date;
$date_file = fopen("./db/latest.updated.date", "r");
if (!feof($date_file)) {
	$line = fgets($date_file);
	$tmparray = preg_split("/\n/", $line);
	$latest_updated_date = $tmparray[0];
}
fclose($date_file);

$date_file = fopen("./db/resolution_summary.txt", "r");
if (!feof($date_file)) {
	$line = fgets($date_file);
	$tmparray = preg_split("/\t/", $line);
	$totol = $tmparray[0];
	if (strlen($totol) > 3) {
		$totol = substr($totol, 0, -3) . "," . substr($totol, -3, 3);
	}
	$n_resolution = $tmparray[1];
	if (strlen($n_resolution) > 3) {
		$n_resolution = substr($n_resolution, 0, -3) . "," . substr($n_resolution, -3, 3);
	}
	$n0 = $tmparray[2];
	if (strlen($n0) > 3) {
		$n0 = substr($n0, 0, -3) . "," . substr($n0, -3, 3);
	}
	$n5 = $tmparray[3];
	if (strlen($n5) > 3) {
		$n5 = substr($n5, 0, -3) . "," . substr($n5, -3, 3);
	}
	$n10 = $tmparray[4];
	if (strlen($n10) > 3) {
		$n10 = substr($n10, 0, -3) . "," . substr($n10, -3, 3);
	}
	$n15 = $tmparray[5];
	if (strlen($n15) > 3) {
		$n15 = substr($n15, 0, -3) . "," . substr($n15, -3, 3);
	}
	$n20 = $tmparray[6];
	if (strlen($n20) > 3) {
		$n20 = substr($n20, 0, -3) . "," . substr($n20, -3, 3);
	}
}
fclose($date_file);


?>

<table width="100%" border="0" cellspacing="0" cellpadding="0">
	<tr>
		<td>
			<table width="700" border="0" align="center" cellpadding="0"
				cellspacing="0">
				<tr valign="top">
					<td width="415">
						<table width="100%" border="0" cellspacing="0" cellpadding="0">
							<tr>
								<td><img src="images/spacer.gif" alt="" width="100" height="35" />
								</td>
							</tr>
							<tr>
								<td>
									<img src="images/surfer.png" alt="" width="83" height="83"
											class="img1" />
									<p class="text6"> <em>EM-SURFER</em> </p>
									<p class="text4"> <a href="submit.php"> Quick Access to Search/Submission Page. </a> </p>
									<p> <strong> A web-based tool for real-time comparison and analysis of Electron Microscopy (EM)
                                                                        density maps.</strong> It compares the shape of EM map isosurfaces, generated using
                                                                        author-recommended contour values. Users can either upload an EM map or choose an existing
                                                                        map in the <a href="http://www.emdatabank.org" target="_blank">EMDataBank</a> and compare it against maps stored
                                                                        in the  EMDataBank. 3D-Zernike Descriptors (3DZD) are utilized for the efficient comparison
                                                                        between EM maps.
                                                                        The web interface will conveniently display a quantitative measure of similarty between
                                                                        the isosurface shape of EM maps, as well as a measure indicating how similar the target
                                                                        isosurface is to the maps retrieved in the search.
                                                                        <em>EM-SURFER</em> provides the users with the option of comparing EM maps using additional
                                                                        contour levels, based on the original values recommended by each map's authors. It also provides
                                                                        filtering features to compare the target only to maps with a similar volume.
									</p>
							   </td>
							</tr>
						</table>
						 <p class="text6"><strong><em>Statistics of latest release  (Last update: <?php echo $latest_updated_date ?>):</em> </strong>
			        </p>
                        <table width="100%" border="0" cellspacing="1" cellpadding="3">
			<!--
		     <tr> 
                        <td class="pdas" bgcolor="#f9f9ff">All  Entries</td>
                        <td bgcolor="#f9f9ff"> 
                        <?php echo "$totol" ?>
                        </td>
                      </tr>-->
                      <tr> 
                        <td class="pdas" bgcolor="#F0F0F0">All  Entries with resolution data</td>
                        <td bgcolor="#F0F0F0"> 
                        <?php echo "$n_resolution" ?>
                        </td>
                      </tr>
                      <tr> 
                        <td class="pdas" bgcolor="#f9f9ff">Entries with resolution <=5&Aring; </td>
                        <td bgcolor="#f9f9ff"> <?php echo $n0 ?></td>
                      </tr>
                      <tr> 
                        <td class="pdas" bgcolor="#F0F0F0">Entries with resolution >5&Aring; and <=10&Aring;</td>
                        <td bgcolor="#F0F0F0"> 
                       <?php echo $n5 ?></td>
                      </tr>
                      <tr> 
                        <td class="pdas" bgcolor="#f9f9ff">Entries with resolution >10&Aring; and <=15&Aring;</td>
                        <td bgcolor="#f9f9ff"> 
                        <?php echo $n10 ?>
                        </td>
                      </tr>
                      <tr>
                        <td class="pdas" bgcolor="#F0F0F0">Entries with resolution >15&Aring; and <=20&Aring;</td> 
                        <td bgcolor="#F0F0F0"><?php echo $n15 ?></td>
                      </tr>
                      <tr>
                        <td class="pdas" bgcolor="#f9f9ff">Entries with resolution >20&Aring;</td> 
                        <td bgcolor="#f9f9ff"> 
                        <?php echo $n20 ?>
                        </td>
                      </tr>
                     
                    </table>
					</td>
					<td width="35">&nbsp;</td>
					
					<td width="250"><p>&nbsp;
							</p>
					  <p class="text6"><strong><em>Latest features:</em> </strong>
				      </p>

						<!-- <marquee behavior="scroll" scrollamount="2" direction="up" > -->
                        <ul>
							<li><strong>Surface representation:</strong> Four isosurface representation options are provided.
							The first one is the isosurface generated at the author-recommended contour level found in the EMDB.
							Three additional options combine the author-recommended contour level with the isosurface created by increasing
							the contour value to represent regions of a higher density, which is usually closer to the core of the
							molecules.</li>
							<br><li><strong>User-uploaded maps:</strong> <em>EM-SURFER</em> now accepts user EM maps as queries for searches.</li>
							<br><li><strong>Volume filter:</strong> Enable this filter to compare only against maps with a volume similar
							to the query map.</li>
							<br><li><strong>Batch mode:</strong> A batch mode is available to query multiple EM maps, for users who would
							like to benchmark their methods against our approach.</li>
						</ul>
                       <!-- </marquee > -->
                       
                    </td>
				</tr>
			</table>
		</td>
	</tr>
</table>

<br>
<br>
<?php
include ('tailer.php');
?>
