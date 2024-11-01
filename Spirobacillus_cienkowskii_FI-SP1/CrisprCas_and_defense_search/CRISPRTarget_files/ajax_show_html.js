function ajaxShowHtml(scriptLink)
{
	//alert(scriptLink);
	//******************************** hide the screen first	
	document.getElementById('submit_button').disabled=true;
	document.getElementById("result_div").innerHTML = "<div style='text-align:center;width:100%;'><img src='../imgs/loading_icon.gif'></img></div>";	
	
				//ajaxCheckLog('/mmdt_cgi-bin/MMDT_AJAX_SCRIPTS/check_log.pl');
	xmlhttpPost_align(scriptLink);
				//ajaxCheckLog('/mmdt_cgi-bin/MMDT_AJAX_SCRIPTS/check_log.pl');
	

	function xmlhttpPost_align(strURL) 
	{
		var xmlHttpReq = false;
		var self = this;
		var i=1;
	   // alert("OK");
		// Mozilla/Safari
		if (window.XMLHttpRequest) 
			{
				self.xmlHttpReq = new XMLHttpRequest();
			}
		// IE
		else if (window.ActiveXObject)
			{
				self.xmlHttpReq = new ActiveXObject("Microsoft.XMLHTTP");
			}
		self.xmlHttpReq.open('POST', strURL, true); //false will make the call synchronous true will make it asynchronous
		self.xmlHttpReq.setRequestHeader('Content-Type', 'application/x-www-form-urlencoded');
		
		self.xmlHttpReq.onreadystatechange = function() 
				{
					if (self.xmlHttpReq.readyState == 4) 
						{

							updatepage(self.xmlHttpReq.responseText);
							//ajaxCheckLog('/mmdt_cgi-bin/MMDT_AJAX_SCRIPTS/check_log.pl');
						
						}
						/* else if (self.xmlHttpReq.readyState == 2) 
								{
							
									//alert("3");
									ajaxCheckLog('/mmdt_cgi-bin/MMDT_AJAX_SCRIPTS/check_log.pl');
											
								}
						*/		 
				}
		self.xmlHttpReq.send(getquerystring());		
	}




	function getquerystring() 
		{
		//	alert("Inside1");
			var form = document.forms["show_html"];
			//alert("Inside2"+form);
		


			
		if(form.spacer_forward_strand.checked)
			{
				var sfs = 1;    	//alert('orf='+orf);
			}
		else{
				var sfs = 0; 
			}	
		if(form.spacer_reverse_strand.checked)
			{
				var srs = 1;    	//alert('orf='+orf);
			}
		else{
				var srs = 0; 
			}
			
								
		
		if(form.score_spacer_match)
			{
				var ssm = form.score_spacer_match.value;    	//alert('orf='+orf);
			}
		if(form.score_spacer_mismatch)
			{
				var ssmismatch = form.score_spacer_mismatch.value;    	//alert('orf='+orf);
			}
			
			
		if(form.flank_length_3p)
			{
				var fl_3p = form.flank_length_3p.value;    	//alert('orf='+orf);
			}
			
		if(form.score_tag_match_3p)
			{
				var stm_3p = form.score_tag_match_3p.value;    	//alert('orf='+orf);
			}
		if(form.score_tag_mismatch_3p)
			{
				var stmismatch_3p = form.score_tag_mismatch_3p.value;    	//alert('orf='+orf);
			}	

		var len=form.selected_PAM_3p.length;    //alert(len);
		var sPAM_3p;
		for (var i=0; i < len; i++)
		   {			  
			  var j=i;			  	  
			   if (form.selected_PAM_3p[i].selected)
				  {
					  
					  sPAM_3p=form.selected_PAM_3p[i].value; 
					  //alert(sPAM);					  
				  }			  
			} //end of for
		
		//alert(sPAM);	
					
		if(form.user_pam_3p)
			{
				var uPAM_3p = form.user_pam_3p.value;    	//alert('orf='+orf);
			}
		else{
				var uPAM_3p = "Not_Specified";
			}
				
		if(form.score_pam_match_3p)
			{
				var spm_3p = form.score_pam_match_3p.value;    	//alert('orf='+orf);
			}
			
		if(form.filter_parameters_3p)
			{
				var fp_3p = form.filter_parameters_3p.value;    	//alert('orf='+orf);
			}
		else{
				var fp_3p =0;
			}		
		if(form.filter_no_mismatch_3p)
			{
				var fnm_3p = form.filter_no_mismatch_3p.value;    	//alert('orf='+orf);
			}
		else{
				var fnm_3p ="NA";
			}		
		if(form.filter_except_base_3p)
			{
				var feb_3p = form.filter_except_base_3p.value;    	//alert('orf='+orf);
			}
		else{
				var feb_3p ="NA";
			}			
			
			


		
		if(form.flank_length_5p)
			{
				var fl_5p = form.flank_length_5p.value;    	//alert('fl_5p='+fl_5p);
			}
			
		if(form.score_tag_match_5p)
			{
				var stm_5p = form.score_tag_match_5p.value;    	//alert('orf='+orf);
			}
		if(form.score_tag_mismatch_5p)
			{
				var stmismatch_5p = form.score_tag_mismatch_5p.value;    	//alert('orf='+orf);
			}				
			
			
		var len2=form.selected_PAM_5p.length;    //alert(len);
		var sPAM_5p;
		for (var i2=0; i2 < len2; i2++)
		   {			  
			  var j2=i2;			  	  
			   if (form.selected_PAM_5p[i2].selected)
				  {
					  
					  sPAM_5p=form.selected_PAM_5p[i2].value; 
					  //alert(sPAM);					  
				  }			  
			} //end of for
		
		//alert(sPAM);	
					
		if(form.user_pam_5p)
			{
				var uPAM_5p = form.user_pam_5p.value;    	//alert('orf='+orf);
			}
		else{
				var uPAM_5p = "Not_Specified";
			}
					
		if(form.score_pam_match_5p)
			{
				var spm_5p = form.score_pam_match_5p.value;    	//alert('orf='+orf);
			}			
			
			
			
		if(form.filter_parameters_5p)
			{
				var fp_5p = form.filter_parameters_5p.value;    	//alert('orf='+orf);
			}
		else{
				var fp_5p =0;
			}		
		if(form.filter_no_mismatch_5p)
			{
				var fnm_5p = form.filter_no_mismatch_5p.value;    	//alert('orf='+orf);
			}
		else{
				var fnm_5p ="NA";
			}		
		if(form.filter_except_base_5p)
			{
				var feb_5p = form.filter_except_base_5p.value;    	//alert('orf='+orf);
			}
		else{
				var feb_5p ="NA";
			}

					
		if(form.group_output.checked)
			{
				var go = 1;    	//alert('orf='+orf);
			}
		else{
				var go = 0;
			}
						
		if(form.cutoff_score)
			{
				var csc = form.cutoff_score.value;    	//alert('csc='+csc);
			}	

		
		
		
		
		if(form.user_spacer_file)
			{
				var upf = form.user_spacer_file.value;    	//alert('orf='+orf);
			}		
			
		if(form.op_intra_features_out)
			{
				var opfo = form.op_intra_features_out.value;    	//alert('orf='+orf);
			}		
		if(form.ref_db_fasta)
			{
				var rdf = form.ref_db_fasta.value;    	//alert('orf='+orf);
			}
		
		if(form.spacer_source_seq_file)
			{
				var sssf = form.spacer_source_seq_file.value;    	//alert('orf='+orf);
			}	
		if(form.id_ref_table)
			{
				var idrt = form.id_ref_table.value;    	//alert('orf='+orf);
			}
			
		if(form.text_report)
			{
				var text_rep = form.text_report.value;    	//alert('orf='+orf);
			}
		if(form.html_report)
			{
				var html_rep = form.html_report.value;    	//alert('html_rep='+html_rep);
			}
		if(form.experiment_title)
			{
				var exp_t = form.experiment_title.value;    //	alert('exp_t='+exp_t);
			}	
		//	+ '&experiment_title=' + exp_t	
		//alert(spm_5p);
						 	
		//if(form.INTERNET_EXPLORER)
		//	{
		//		var bName = navigator.appName;
		//		//var bVer = integer(navigator.appVersion);			   
		//		if (bName == "Microsoft Internet Explorer"){var ieB=1;}
		//		else{var ieB=0;}
				
				//var ieB = form.INTERNET_EXPLORER.value;	alert('ieB='+ieB);
		//	}		
				
		var qstr="";
		
		qstr =	'score_spacer_match=' + ssm + '&user_spacer_file=' + upf + '&text_report=' + text_rep + '&html_report=' + html_rep + '&spacer_forward_strand=' + sfs + '&spacer_reverse_strand=' + srs + '&cutoff_score=' + csc +'&score_spacer_mismatch=' + ssmismatch + '&flank_length_3p=' + fl_3p + '&score_tag_match_3p=' + stm_3p + '&score_tag_mismatch_3p=' + stmismatch_3p + '&' + 'selected_PAM_3p=' + sPAM_3p + '&' + 'user_pam_3p=' + uPAM_3p + '&' + 'score_pam_match_3p=' + spm_3p + '&filter_parameters_3p='+ fp_3p + '&filter_no_mismatch_3p=' + fnm_3p + '&filter_except_base_3p=' + feb_3p + '&flank_length_5p=' + fl_5p + '&score_tag_match_5p=' + stm_5p + '&score_tag_mismatch_5p=' + stmismatch_5p + '&' + 'selected_PAM_5p=' + sPAM_5p + '&' + 'user_pam_5p=' + uPAM_5p + '&'+ 'score_pam_match_5p=' + spm_5p +  '&group_output=' + go + '&op_intra_features_out=' + opfo + '&ref_db_fasta=' + rdf + '&spacer_source_seq_file=' + sssf + '&experiment_title=' + exp_t + '&id_ref_table=' + idrt;
		
		//qstr ='ORIGINAL_REF_FILE=' + orf + '&' + 'LAST_PROCESSED_FILE=' + lpf + '&' + 'HIGHEST_FILE_VERSION_NO=' + hfv + '&' + 'LOG_FILE=' + lf + '&' + 'INTERNET_EXPLORER=' + ieB;  // NOTE: no '?' before querystring
	   
		//alert(qstr);
	 
		return qstr;
		}

	function updatepage(str)
		{
			//var ih=document.getElementById("LPF").innerHTML;
			
			//alert(ih);
			/* if(str.match(/PRINT_START/g))
				{
					//alert(str);
					//document.getElementById("logs").style.display ="none";
					//document.getElementById("result_div").innerHTML = "";
					document.getElementById("result_div").innerHTML = str;
					//document.getElementById("html_content").innerHTML = str;
				}
			else{	
					alert(str);
					document.getElementById("logs_content").innerHTML = str;
				}
			*/
			
			document.getElementById("result_div").innerHTML = "";	
			document.getElementById("result_div").innerHTML = str;	
			document.getElementById('submit_button').disabled=false;
			//document.getElementById('align_button').innerHTML="Align";
			//document.getElementById("align_loader").style.display ="none";			
			//document.getElementById('BlockDiv').style.display='none';
			
			//ajaxCheckLog('/mmdt_cgi-bin/MMDT_AJAX_SCRIPTS/check_log.pl');
		}

}
