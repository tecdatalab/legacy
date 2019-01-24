var xmlHttp;
var rsmdArray = new Array();

// HTML ID of all VisGrid buttons used
//var visGridButtons = ["cavitytop", "protrusiontop", "flattop", "cleartop","cavitybottom", "protrusionbottom", "flatbottom", "clearbottom"];
var visGridButtons = ["cavitytop", "protrusiontop", "flattop", "pockettop","cavitybottom", "protrusionbottom", "flatbottom", "pocketbottom"];

/*
 * This function returns the number of elements that are being displayed
 * based on the value of the selection list called "num"
 */
function howManyAreDisplayed()
{
	var selectList = document.getElementsByName("num");
	return selectList[0].value;
}

function calculateRMSD(query,j)
{
	var i;
	xmlHttp = getXmlHttpObject();
	if (xmlHttp == null)
	{
		alert("AJAX not supported");
	}
	var url = "calculateRMSD.cgi";
	url = url + "?q=" + query;
	for (i = 1; i <= howManyAreDisplayed(); i++)
	{
		url = url + "&rmsd" + i + "=";
		var rmsdCheckbox = document.getElementById("rmsdbox" + i);
		// if the element wasn't found then it has already been replaced by a button
		// skip this entry
		if(rmsdCheckbox == null)
		{
			continue;
		}
		//if (eval("document.display.rmsdbox" + i + ".checked") == true)
		if (i==j)
		{
			url = url + rmsdCheckbox.value;
			rmsdCheckbox.checked = false;
			// replace the checkbox by a button to display the RMSD difference in JMol
			document.getElementById("rmsd" + i).innerHTML = "<img src=\"../loadingbar2.gif\" vspace=\"2\" hspace=\"2\" />";
		}
		else 
		{
			url = url + "0";
		}
	}
	// before starting the operation, disable functions
	enableFunctions(false);
	xmlHttp.onreadystatechange = stateChanged;
	xmlHttp.open("GET", url, true);
	xmlHttp.send(null);
}

function calculateRMSDCustom(query,j)
{
	var i;
	xmlHttp = getXmlHttpObject();
	if (xmlHttp == null)
	{
		alert("AJAX not supported");
	}
	var url = "calculateRMSDCustom.cgi";
	url = url + "?q=" + query;
	for (i = 1; i <= howManyAreDisplayed(); i++)
	{
		url = url + "&rmsd" + i + "=";
		var rmsdCheckbox = document.getElementById("rmsdbox" + i);
		// if the element wasn't found then it has already been replaced by a button
		// skip this entry
		if(rmsdCheckbox == null)
		{
			continue;
		}
		
		//if (eval("document.display.rmsdbox" + i + ".checked") == true)
		if (i==j)
		{
			url = url + rmsdCheckbox.value;
			rmsdCheckbox.checked = false;
			// replace the checkbox by a button to display the RMSD difference in JMol
			document.getElementById("rmsd" + i).innerHTML = "<img src=\"../loadingbar2.gif\" vspace=\"2\" hspace=\"2\" />";
		}
		else
		{
			url = url + "0";
		}
	}
	// before starting the operation, disable functions
	enableFunctions(false);
	xmlHttp.onreadystatechange = stateChanged;
	xmlHttp.open("GET", url, true);
	xmlHttp.send(null);
}


function visgrid(query, type)
{
	//document.jmol.script('select all; color white');
	var resultTextArea = document.getElementById("statusbar");
	if (vginit == 1)
	{
		if (type == "cavity")
		{
			//resultJmol = jmolScriptWait(vgcavity);
			document.jmol.script(vgcavity);
			resultTextArea.value = cavRes;
		}
		if (type == "protrusion")
		{
			//resultJmol = jmolScriptWait(vgprotrusion);
			document.jmol.script(vgprotrusion);
			resultTextArea.value = protRes;
		}
		if (type == "flat")
		{
			//resultJmol = jmolScriptWait(vgflat);
			document.jmol.script(vgflat);
			resultTextArea.value = flatRes;
		}
		vginit = 0;
		return;
	}
	xmlHttp = getXmlHttpObject();
	if (xmlHttp == null)
	{
		alert("AJAX not supported");
	}
	document.getElementById("visgridload1").innerHTML = "<img src=\"../loadingbar2.gif\" vspace=\"2\" hspace=\"2\" />";
	document.getElementById("visgridload2").innerHTML = "<img src=\"../loadingbar2.gif\" vspace=\"2\" hspace=\"2\" />";
	var url = "jmol_visgrid.cgi?pdbid=" + query + "&type=" + type;
	enableFunctions(false);
	xmlHttp.onreadystatechange = stateChanged;
	xmlHttp.open("GET", url, true);
	xmlHttp.send(null);
}


function visgridCustom(query, type)
{
	//document.jmol.script('select all; color white');
	var resultTextArea = document.getElementById("statusbar");
	if (vginit == 1)
	{
		if (type == "cavity")
		{
			//resultJmol = jmolScriptWait(vgcavity);
			document.jmol.script(vgcavity);
			resultTextArea.value = cavRes;
		}
		if (type == "protrusion")
		{
			//resultJmol = jmolScriptWait(vgprotrusion);
			document.jmol.script(vgprotrusion);
			resultTextArea.value = protRes;
		}
		if (type == "flat")
		{
			//resultJmol = jmolScriptWait(vgflat);
			document.jmol.script(vgflat);
			resultTextArea.value = flatRes;
		}
		vginit = 0;
		return;
	}
	xmlHttp = getXmlHttpObject();
	if (xmlHttp == null)
	{
		alert("AJAX not supported");
	}
	document.getElementById("visgridload1").innerHTML = "<img src=\"../loadingbar2.gif\" vspace=\"2\" hspace=\"2\" />";
	document.getElementById("visgridload2").innerHTML = "<img src=\"../loadingbar2.gif\" vspace=\"2\" hspace=\"2\" />";
	var url = "jmol_visgridCustom.cgi?pdbid=" + query + "&type=" + type;
	enableFunctions(false);
	xmlHttp.onreadystatechange = stateChanged;
	xmlHttp.open("GET", url, true);
	xmlHttp.send(null);
}

function ligsitecs(query, type)
{
	var resultTextArea = document.getElementById("statusbar");
	if (vginit == 1)
	{
		if (type == "pocket")
		{
			document.jmol.script(vgcavity);
			resultTextArea.value = cavRes;
		}
		vginit = 0;
		return;
	}
	xmlHttp = getXmlHttpObject();
	if (xmlHttp == null)
	{
		alert("AJAX not supported");
	}
	document.getElementById("ligsitecsload1").innerHTML = "<img src=\"../loadingbar2.gif\" vspace=\"2\" hspace=\"2\" />";
	document.getElementById("ligsitecsload2").innerHTML = "<img src=\"../loadingbar2.gif\" vspace=\"2\" hspace=\"2\" />";
	var url = "jmol_ligsitecs.cgi?pdbid=" + query + "&type=" + type;
	enableFunctions(false);
	xmlHttp.onreadystatechange = stateChanged;
	xmlHttp.open("GET", url, true);
	xmlHttp.send(null);
}

function ligsitecsCustom(query, type)
{
	var resultTextArea = document.getElementById("statusbar");
	if (vginit == 1)
	{
		if (type == "pocket")
		{
			document.jmol.script(vgcavity);
			resultTextArea.value = cavRes;
		}
		vginit = 0;
		return;
	}
	xmlHttp = getXmlHttpObject();
	if (xmlHttp == null)
	{
		alert("AJAX not supported");
	}
	document.getElementById("ligsitecsload1").innerHTML = "<img src=\"../loadingbar2.gif\" vspace=\"2\" hspace=\"2\" />";
	document.getElementById("ligsitecsload2").innerHTML = "<img src=\"../loadingbar2.gif\" vspace=\"2\" hspace=\"2\" />";
	var url = "jmol_ligsitecsCustom.cgi?pdbid=" + query + "&type=" + type;
	enableFunctions(false);
	xmlHttp.onreadystatechange = stateChanged;
	xmlHttp.open("GET", url, true);
	xmlHttp.send(null);
}

function jmolCEPose(i,query_id,result_id,rmsd, coverage)
{

	// Run CEPose 
	xmlHttp = getXmlHttpObject();
	if (xmlHttp == null)
	{
		alert("AJAX not supported");
	}
	var url = "jmolCEPose.cgi";
	i++;
	url = url + "?n=" + i + "&q=" + query_id + "&r=" + result_id + "&rmsd=" + rmsd + "&coverage=" + coverage;
	document.getElementById("rmsd" + i).innerHTML = "<img src=\"../loadingbar2.gif\" vspace=\"2\" hspace=\"2\" />";
	
	// before starting the operation, disable RMSD functions
	enableFunctions(false);

	xmlHttp.onreadystatechange = stateChanged;
	xmlHttp.open("GET", url, true);
	xmlHttp.send(null);
}

// When invoking RMSD calculation or display operations all other RMSD functions should be disabled. After the operation
// finishes, they should be re-established. This function enables or disables all rmsd checkboxes and buttons in these situations
function enableRMSDFunctions(isEnabled)
{
	var inconsistency_found = false;
	var inconsistency_id;
	for (i = 1; i <= howManyAreDisplayed(); i++)
	{
		var rmsdCheckbox = document.getElementById("rmsdbox" + i);
		// the button collection should contain either 0 or 1 buttons
		var rmsdButtonCollection = document.getElementsByName("rmsdbutton" + i);
		// if neither a checkbox nor button is found, then ignore it (although this shouldn't happen)
		if(rmsdCheckbox != null) //disable a checkbox
		{
			rmsdCheckbox.disabled = !isEnabled;
		}
		else if(rmsdButtonCollection.length > 0) // we have a button: disable it
		{
			rmsdButtonCollection[0].disabled = !isEnabled;
		}
		// else if neither a checkbox nor a button is found, then the calculatiion failed (although this shouldn't happen)
		else
		{
			if(isEnabled) // if we are re-enabling functions then log the inconsistency
			{
				inconsistency_found = true;
				inconsistency_id = i;
			}
		}
	}
	if(inconsistency_found) // avoid loading bar loops by replacing inconsistencies by N/A
	{
		document.getElementById("rmsd" + inconsistency_id).innerHTML = "Rmsd: N/A";
	}

}

// Similar to the de-activation of RMSD functions, when an operation is being processed
// the VisGrid functions will be deactivated to avoid problems because of parallel processing.
// If a false is passed as parameter, then all of the VisGrid buttons will be disabled, while a true
// will re-enable them
function enableVisGridFunctions(isEnabled)
{
	for(buttonIndex=0; buttonIndex < visGridButtons.length; buttonIndex++)
	{
		var button = visGridButtons[buttonIndex];
		document.getElementById(button).disabled = !isEnabled;
	}
}

// Wrapper method for other more specific functions that enable/disable the different parts of the interface
// At this moment it only takes into account RMSD and VisGrid, but if we add more functionality later on
// a new function to enable/disable that section should be created, and it should be invoked from here
function enableFunctions(isEnabled)
{
	enableRMSDFunctions(isEnabled);
	enableVisGridFunctions(isEnabled);
}

function stateChanged()
{
	if (xmlHttp.readyState == 4)
	{
//		document.getElementById("testA").innerHTML = xmlHttp.responseText;
		var ajaxArray = xmlHttp.responseText.split("`;");
		if (ajaxArray[0] == "rmsd")
		{
			var query_id = ajaxArray[51];
			var result_id = ajaxArray[52];
                        var coverage = ajaxArray[53];
			var result_found = false;
			for (i = 1; i <= 25; i++)
			{
				if (ajaxArray[i] != 0)
				{

                                        if (ajaxArray[i] == "N/A") {
                        			ce_rmsd_str = ajaxArray[i];
                        		} else {
						ce_rmsd_str = ajaxArray[i] + "Å";
                        		}
					document.getElementById("rmsd" + i).innerHTML = "<input type=\"button\" name=\"rmsdbutton" + i +
											"\" value=\"Rmsd\" onClick=\"javascript:jmolCEPose(\'" +
											(i - 1) + "\',\'" + query_id + "\',\'" + result_id + "\',\'" +
											ajaxArray[i] + "\',\'" + coverage + "\')\">: " +
											"<font color=\"blue\"><u>" +
											"<a href=\"tmp/ce_" + result_id + "_" + i + ".txt\" target=\"_blank\">" +
											ce_rmsd_str + coverage +
											"</a></u></font>";
				}
			}
			
			for (i = 26; i <= 50; i++)
			{
				if (ajaxArray[i] != 0)
				{
					ceArray[i - 26] = ajaxArray[i];
//					alert(i + ceArray[i - 26]);
				}
			}
			// re-enable checkboxes and buttons, because the operation finished
			enableFunctions(true);
		}
		if (ajaxArray[0] == "visgrid")
		{
			document.getElementById("visgridload1").innerHTML = "";
			document.getElementById("visgridload2").innerHTML = "";
//			document.jmol.script(ajaxArray[1]);
			vginit = 1;
			vgcavity = ajaxArray[1];
			vgprotrusion = ajaxArray[2];
			vgflat = ajaxArray[3];
			cavRes = ajaxArray[4];
			protRes = ajaxArray[5];
			flatRes = ajaxArray[6];
			visgrid("", ajaxArray[7]);
			enableFunctions(true);
		}

 //add 
                if (ajaxArray[0] == "ligsitecs")
		{
			document.getElementById("ligsitecsload1").innerHTML = "";
			document.getElementById("ligsitecsload2").innerHTML = "";
// document.jmol.script(ajaxArray[1]);
			vginit = 1;
			vgcavity = ajaxArray[1];
			cavRes = ajaxArray[2];
			ligsitecs("", ajaxArray[3]);
			enableFunctions(true);
		}

		if (ajaxArray[0] == "jmolcepose") 
		{
			ce_query = ajaxArray[1]
			ce_result = ajaxArray[2];
			ce_i = ajaxArray[3];
			ce_rmsd = ajaxArray[4];
			if (ce_rmsd == "N/A") {
                        	ce_rmsd_str = ce_rmsd;
                        } else {
				ce_rmsd_str = ce_rmsd + "Å";
                        }
                        coverage = ajaxArray[5];

		
			// Update Jmol
			jmcmd = "zap; load ../cgi-bin/tmp/cepose\_" + ce_result + "\_" + ce_i + ".pdb" +  "; select all; spacefill off; wireframe off; backbone 100; select :A; color red; select :B; color blue;";
			//resultJmol = jmolScriptWait(jmcmd);
			document.jmol.script(jmcmd);
			
			document.getElementById("rmsd" + ce_i).innerHTML = "<input type=\"button\" name=\"rmsdbutton" + ce_i +
									"\" value=\"Rmsd\" onClick=\"javascript:jmolCEPose(\'" + (ce_i - 1) +
									"\',\'" + ce_query + "\',\'" + ce_result + "\',\'" + ce_rmsd + "\')\">: " + 
									"<font color=\"blue\"><u>" +
									"<a href=\"tmp/ce_" + ce_result + "_" + ce_i + ".txt\" target=\"_blank\">" +
									ce_rmsd_str + coverage +
									"</a></u></font>";
									
			// re-enable checkboxes and buttons, because the operation finished
			enableFunctions(true);
		}
	}
}

function getXmlHttpObject()
{
	var xmlHttp=null;
	try
	{
		xmlHttp = new XMLHttpRequest();
	}
	catch (e)
	{
		try
		{
			xmlHttp = new ActiveXObject("Msxml2.XMLHTTP");
		}
		catch (e)
		{
			xmlHttp = new ActiveXObject("Microsoft.XMLHTTP");
		}
	}
	return xmlHttp;
}
